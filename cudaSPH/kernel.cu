
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <stdlib.h>
#include <crtdbg.h>
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>
#include <algorithm>

//#include "particles.h"
#include "shader.h"

#include <vector>
#include <math.h>
//#include <vld.h>

#define _CRTDBG_MAP_ALLOC

#define kParticlesCount 1024
#define kBorderSideParticlesCount 300
#define kBorderParticlesCount (4*kBorderSideParticlesCount)
#define kWindowWidth 640
#define kWindowHeight 480
#define kPi 3.14159265359
#define g -9.81
#define kFrameRate 60
#define kSubSteps 7

#define kOffset 0.90f

//#define kDt ((1.0f/kFrameRate) / kSubSteps)

#define averageParticles 20
//#define interactionRadius sqrt(averageParticles/(kParticlesCount*kPi))
#define interactionRadius 0.1f
#define IR2 interactionRadius*interactionRadius
#define cellSize (2.0f*interactionRadius)

void advance();
void render();
void glInit();
void particlesInit();
void borderParticlesInit();
void drawGrid();
void updateGrid();
void updatNeighbours();
void createDrawablePoints();
void calculatePressure();
void calulateForces();
float calculateMass();


#define kWdeafult (315/(64*kPi*pow(interactionRadius,9))); 
#define kWgradPressure -(45/(kPi*pow(interactionRadius,6)));
#define kWlaplacianViscosity (45/(kPi*pow(interactionRadius,6)));
#define kWgradDefult -(945/(32*kPi*pow(interactionRadius,9)));
#define kWlaplacianDefult -(945/(32*kPi*pow(interactionRadius,9)));

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

const float kViewScale =  2.0f;
//const float interactionRadius =  0.05f;
//const float cellSize = 2*interactionRadius;



#define kDt  0.0005f;
const int kCellCount = 100;
const float restDensity = 988.0f;
const int kstiffnes = 20;
const float surfaceTension = 0.0728f;
const float viscosityConstant = 3.5f;

float particleMass;
float surfaceLimit = sqrt(restDensity/averageParticles);
float accelerationX;
float accelerationY;

//Every force
float   pressureForcex, pressureForcey, viscosityForcex,
        viscosityForcey,normalx, normaly, gradNormal,
        surfaceTensionForcex,surfaceTensionForcey, gravity;

float dx, dy, distance2;


struct particle
{
    float m_x;
    float m_y;
    float m_u; //x-velocity
    float m_v; //y-velocity

    float m_massDensity;
    float m_pressure;

    float m_mass;

    particle* next;
};

#define kMaxNeighbourCount 64
struct Neighbours
{
	//thrust::device_vector<particle*> particlesVec[kMaxNeighbourCount];
	const particle* particles[kMaxNeighbourCount];
    float r2[kMaxNeighbourCount];
    size_t count;
};

__shared__ Neighbours self_label;

struct point
{
    float x;
    float y;
};


particle particles[kParticlesCount];

particle borderParticles[kBorderParticlesCount];
const size_t kGridWidth = (size_t)(2.0 / cellSize);
const size_t kGridHeight = (size_t)(2.0 / cellSize);

const size_t kGridCellCount = kGridWidth * kGridHeight;
size_t gridCoords[kParticlesCount*2];
particle* gridArr[kGridCellCount];
std::vector<particle*> grid;

std::vector<point> drawablePoints;
std::vector<point> DEBUG_CORNER;

point acceleration[kParticlesCount];
point prevAcceleration[kParticlesCount];
float vhx[kParticlesCount];
float vhy[kParticlesCount];
bool firstIteration = true;

Neighbours neighbours[kParticlesCount];
GLuint programID = 0;

//Global CUDA arrays
__device__ particle  d_particles[kParticlesCount];
__device__ particle* d_grid[kGridCellCount];
thrust::device_vector<particle*> ddgrid;
 //thrust::device_vector<particle*> d_grid;
__device__ size_t d_gridCoords[2*kParticlesCount];
__device__ Neighbours d_neighbours[kParticlesCount];
__device__ particle d_borderParticles[kBorderParticlesCount];
__device__ float d_vhx[kParticlesCount];
__device__ float d_vhy[kParticlesCount]; 
__device__ bool d_firstIteration = true;

__global__ void updateGridDevice()
{
    const size_t d_kGridWidth = (size_t)(2.0 / cellSize);
    const size_t d_kGridHeight = (size_t)(2.0 / cellSize);
    const size_t d_kGridCellCount = d_kGridWidth * d_kGridHeight;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    
        particle& pi = d_particles[i];

        int x = (1 + pi.m_x)/cellSize;
        int y = (1 + pi.m_y)/cellSize;


        if (x < 1)
            x = 1;
        else if (x > d_kGridWidth-2)
            x = d_kGridWidth-2;

        if (y < 1)
            y = 1;
        else if (y > d_kGridHeight-2)
            y = d_kGridHeight-2;

        pi.next = d_grid[x+y*d_kGridWidth];
        d_grid[x+y*d_kGridWidth] = &pi;

        d_gridCoords[i*2] = x;
        d_gridCoords[i*2+1] = y;
}

__global__ void updatNeighboursDevice()
{
	const size_t d_kGridWidth = (size_t)(2.0 / cellSize);
	int i = blockIdx.x * blockDim.x + threadIdx.x;

        particle& pi = d_particles[i];

        int x = (1 + pi.m_x)/cellSize;
        int y = (1 + pi.m_y)/cellSize;

        size_t gi = d_gridCoords[i*2];
        size_t gj = d_gridCoords[i*2+1]*d_kGridWidth;
        
        d_neighbours[i].count = 0;
        
        //Loop over border
        for(size_t j = 0; j < kBorderParticlesCount; j++)
        {
            particle bp = d_borderParticles[j];
            float pm = bp.m_mass;

            float dx = pi.m_x - bp.m_x;
            float dy = pi.m_y - bp.m_y;
            float distance2 = dx*dx + dy*dy;

            if(distance2 < IR2)
            {
                if(d_neighbours[i].count < kMaxNeighbourCount)
                {
                    d_neighbours[i].particles[d_neighbours[i].count] = &bp;
                    d_neighbours[i].r2[d_neighbours[i].count] = distance2;
                    ++d_neighbours[i].count;
                }
            }
        }
        //loop over cells
        for (int ni=gi-1; ni<=gi+1; ++ni)
        {
            for (int nj=gj-d_kGridWidth; nj<=gj+d_kGridWidth; nj+=d_kGridWidth)
            {
                //loop over neighbors
                for (particle* ppj=d_grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
                {
                    //do fancy math
                    //std::cout << "ppj x: " << ppj->m_x << std::endl;
                    float dx = pi.m_x - ppj->m_x;
                    float dy = pi.m_y - ppj->m_y;
                    float distance2 = dx*dx + dy*dy;

                    if(distance2 < IR2)
                    {
                        //Density
                        //massDensity += particleMass*kWdeafult* (IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2);

                        if(d_neighbours[i].count < kMaxNeighbourCount)
                        {
                            d_neighbours[i].particles[d_neighbours[i].count] = ppj;
                            d_neighbours[i].r2[d_neighbours[i].count] =distance2;
                            ++d_neighbours[i].count;
                        }
                    }
                }
            }
        }
}

__global__ void updateGridDevice2()
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	particle& pi = d_particles[i];

	const size_t d_kGridWidth = (size_t)(2.0 / cellSize);
    const size_t d_kGridHeight = (size_t)(2.0 / cellSize);
    const size_t d_kGridCellCount = d_kGridWidth * d_kGridHeight;

	int x = (1 + pi.m_x)/cellSize;
    int y = (1 + pi.m_y)/cellSize;

	if (x < 1)
            x = 1;
        else if (x > d_kGridWidth-2)
            x = d_kGridWidth-2;

	if (y < 1)
            y = 1;
        else if (y > d_kGridHeight-2)
            y = d_kGridHeight-2;

	pi.next = d_grid[x+y*d_kGridWidth];
    d_grid[x+y*d_kGridWidth] = &pi;

    d_gridCoords[i*2] = x;
    d_gridCoords[i*2+1] = y;

	//pi.m_u = d_gridCoords[i*2];
	//pi.m_v = d_gridCoords[i*2+1];

	pi.m_u = d_grid[x+y*d_kGridWidth]->m_mass;
	//pi.m_v = pi.next->m_mass;
}

__global__ void updateNeighboursDevice(){
	const size_t d_kGridWidth = (size_t)(2.0 / cellSize);
	int i = blockIdx.x * blockDim.x + threadIdx.x;

    particle& pi = d_particles[i];

	int x = (1 + pi.m_x)/cellSize;
	int y = (1 + pi.m_y)/cellSize;

	size_t gi = d_gridCoords[i*2];
	size_t gj = d_gridCoords[i*2]*d_kGridWidth;
	/*
	for(int ni = gi-1; ni<=gi+1; ++ni)
	{
		for(int nj = gj-d_kGridWidth; nj<=gj+d_kGridWidth; nj+=d_kGridWidth)
		{
			for(particle* ppj=d_grid[ni+nj]; NULL !=ppj; ppj=ppj->next)
			{

			}
		}
	}
	*/
	//particle* pj = d_grid[i];  
	pi.m_u = d_grid[i]->m_mass;
	//pi.m_v = 0.0f;
	//pi.m_v = gj;
}
__global__ void calculatePressureDevice(){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	particle& pi = d_particles[i];
	int count1 = 0;
	float massDensity = 0.0f;
	int count2 = 0;
	for(int j=0; j< kParticlesCount; j++){
		
		particle& pj = d_particles[j];
		float dx = pi.m_x - pj.m_x;
		float dy = pi.m_y - pj.m_y;
		float distance2 = dx*dx + dy*dy;
		float mass = pj.m_mass;
		
		if(distance2 < IR2){
			count1++;
			massDensity += mass*(IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2)*kWdeafult;
			if(distance2 != 0){
				count2++;		
			}
		}
	}
	
	pi.m_massDensity = massDensity;
	pi.m_pressure = kstiffnes*(massDensity - 988);
	//pi.m_u = count1;
	//pi.m_v = count2;
}

__global__ void calclateForceDevice(){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	particle& pi = d_particles[i];
	
	float pressureForcex = 0.0f; 
    float pressureForcey = 0.0f;
    float viscosityForcex = 0.0f;
    float viscosityForcey = 0.0f;
    float normalx = 0.0f;
    float normaly = 0.0f;
    float gradNormal = 0.0f;
    float surfaceTensionForcex = 0.0f;
    float surfaceTensionForcey = 0.0f;

	float mdi = pi.m_massDensity;

	int count1 = 0;
	int count2 = 0;
	float k = 0.0f;
	for(int j=0; j< kParticlesCount; j++){
		
		particle& pj = d_particles[j];
		float dx = pi.m_x - pj.m_x;
		float dy = pi.m_y - pj.m_y;
		float distance2 = dx*dx + dy*dy;
		float mass = pj.m_mass;
		
		if(distance2 < IR2){
			float mdj = pj.m_massDensity;

			normalx += (mass/mdj)*(IR2-distance2)*(IR2-distance2)*dx*kWgradDefult;
            normaly += (mass/mdj)*(IR2-distance2)*(IR2-distance2)*dy*kWgradDefult;
			gradNormal += (mass/mdj)*(IR2-distance2)*(IR2-distance2)*(3*IR2-7*distance2)*kWlaplacianDefult;
			
			count1++;
			//massDensity += mass*(IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2)*kWdeafult;
			if(distance2 != 0){
				float distance = sqrt(distance2);
				float velocityDiffu = pj.m_u - pi.m_u;
                float velocityDiffv = pj.m_v - pi.m_v;

				float pressi = pi.m_pressure;
				float pressj = pj.m_pressure;

				pressureForcex += ((pressi/(mdi*mdi))+(pressj/(mdj*mdj)))*mass*(interactionRadius-distance)*(interactionRadius-distance)*(dx/distance)*kWgradPressure;
                pressureForcey += ((pressi/(mdi*mdi))+(pressj/(mdj*mdj)))*mass*(interactionRadius-distance)*(interactionRadius-distance)*(dy/distance)*kWgradPressure;

				viscosityForcex += velocityDiffu * (mass/mdj) * (interactionRadius-distance)*kWlaplacianViscosity;
                viscosityForcey += velocityDiffv * (mass/mdj) * (interactionRadius-distance)*kWlaplacianViscosity;

			}
		}
	}

	float normalLenght = 1/sqrt(normalx*normalx + normaly*normaly);
    if(normalLenght > 7){
		surfaceTensionForcex = - 0.0728f  * gradNormal * normalx * normalLenght;
        surfaceTensionForcey = - 0.0728f  * gradNormal * normaly *normalLenght ;
    }

	pressureForcex = -mdi*pressureForcex;
	pressureForcey = -mdi*pressureForcey;

	viscosityForcex = 3.5f * viscosityForcex;
	viscosityForcey = 3.5 * viscosityForcey;

	float accX = (pressureForcex + viscosityForcex + surfaceTensionForcex)/mdi;
	float accY = ((pressureForcey + viscosityForcey + surfaceTensionForcey)/mdi);
	
	/*
	if(d_firstIteration){
		d_vhx[i] = pi.m_u + 0.5*accX*kDt;
        d_vhy[i] = pi.m_v + 0.5*accY*kDt;
                        
        pi.m_u += accX*kDt;
        pi.m_v += accY*kDt;
		
		pi.m_x += d_vhx[i]*kDt;
        pi.m_y += d_vhy[i]*kDt;

        d_firstIteration = false;
	}else{
		d_vhx[i] += accX*kDt;
        d_vhy[i] += accY*kDt;
		pi.m_u = d_vhx[i] + 0.5*accX*kDt;
		pi.m_v = d_vhy[i] + 0.5*accY*kDt;

		pi.m_x += d_vhx[i]*kDt;
		pi.m_y += d_vhy[i]*kDt;
	}
	*/
}


int main () {
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit ()) {
    fprintf (stderr, "ERROR: could not start GLFW3\n");
    return 1;
  } 


  GLFWwindow* window = glfwCreateWindow (kWindowWidth, kWindowHeight, "Here be fluids", NULL, NULL);
  if (!window) 
  {
    fprintf (stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent (window);
                                  
  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  glewInit ();

  // get version info
  const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
  const GLubyte* version = glGetString (GL_VERSION); // version as a string
  printf ("Renderer: %s\n", renderer);
  printf ("OpenGL version supported %s\n", version);
  

  grid.resize(kGridCellCount);
  particlesInit();
  borderParticlesInit();
  updateGrid();
  drawablePoints.resize(kParticlesCount + kBorderParticlesCount);


  particleMass = calculateMass();
  std::cout << "Mass: "<< particleMass << std::endl;
  
  cudaMemcpyToSymbol(d_particles, particles, kParticlesCount*sizeof(particle));
  cudaMemcpyToSymbol(borderParticles, borderParticles, kBorderParticlesCount*sizeof(particle));
  //updateGridDevice<<< 1,kParticlesCount >>>();

  glInit();

  std::cout << "interaction radius: " << interactionRadius << std::endl;
  std::cout << "cellSize: " << cellSize << std::endl;
  std::cout << "grid cell count: " << kGridCellCount << std::endl;
  std::cout << "kgridwidth: " << kGridWidth << std::endl;
  //std::cout << "kDt: " << kDt << std::endl;


  double t = 0.0;
  double currentTime = glfwGetTime();
  double accumulator = 0.0;

  while (!glfwWindowShouldClose (window)) 
  {

	  /*for(int i = 0; i < kSubSteps; ++i)
      {
		updateGridDevice <<< 1, kParticlesCount >>>();
		updatNeighboursDevice <<< 1,kBorderParticlesCount >>> ();
		calculatePressureDevice<<< 16,64 >>>();
		calclateForceDevice<<<16, 64 >>>();
      }*/

	 updateGridDevice2<<< 1,kParticlesCount >>>();
	 //updateNeighboursDevice <<< 1,kParticlesCount>>>();

	  cudaMemcpyFromSymbol(particles, d_particles, kParticlesCount*sizeof(particle));

		for(int i=0; i<kParticlesCount; i++)
		{
			std::cout << "count1 = " << particles[i].m_u << std::endl;
			std::cout << "count2 = " << particles[i].m_v << std::endl;
		}
		
		
      createDrawablePoints();
      render();
	 

      glfwPollEvents ();
      // put the stuff we've been drawing onto the display
      glfwSwapBuffers (window);
  }


  // close GL context and any other GLFW resources
  glfwTerminate();
  return 0;
}

void glInit()
{
     programID = LoadShader( "default.vert", "flat.frag" );
    


    point topleft;
    point topright;
    point bottomleft;
    point bottomright;
    
    topleft.x = -1.0f;
    topleft.y = 1.0f;

    topright.x = 1.0f;
    topright.y = 1.0f;

    bottomleft.x = -1.0f;
    bottomleft.y = -1.0f;

    bottomright.x = 1.0f;
    bottomright.y = -1.0f;

    DEBUG_CORNER.push_back(topleft);
    DEBUG_CORNER.push_back(topright);
    DEBUG_CORNER.push_back(bottomleft);
    DEBUG_CORNER.push_back(bottomright);

    vbo = 0;
    glGenBuffers (1, &vbo);
    

    vao = 0;
    glGenVertexArrays (1, &vao);
    
    


}
void createDrawablePoints()
{
	
	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < kParticlesCount; ++i)
	{
		point p;
		p.x = particles[i].m_x;
		p.y = particles[i].m_y;
		drawablePoints[i] = p;
	}
	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < kBorderParticlesCount; ++i)
	{
		point p;
		p.x = borderParticles[i].m_x;
		p.y = borderParticles[i].m_y;
		drawablePoints[kParticlesCount + i] = p;
	}
}

void render()
{

    glClearColor(0.05f, 0.05f, 0.05f, 1);
	glClear (GL_COLOR_BUFFER_BIT);
	glUseProgram (programID);
    
	glBindBuffer (GL_ARRAY_BUFFER, vbo);
	glBufferData (GL_ARRAY_BUFFER, (kParticlesCount + kBorderParticlesCount) * sizeof(point), &drawablePoints[0], GL_STATIC_DRAW);
	glBindVertexArray (vao);
    //glBindVertexArray (vao);
	glEnableVertexAttribArray (0);
	glBindBuffer (GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer (0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    
 
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 640, 0, 480, 0, 1);
 
	//Draw points as smooth balls (with AA)
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glPointSize(5.0f);

    // draw points from the currently bound VAO with current in-use shader
    glDrawArrays (GL_POINTS, 0, kParticlesCount + kBorderParticlesCount);
	//glDisableVertexAttribArray(0);
}

void particlesInit()
{

	//std::mt19937 eng((std::random_device())());
	//std::uniform_real_distribution<> pos_dist(-0.001,0.001);

	int rowcolSize = sqrt(kParticlesCount);

  for(int particleIndexRow = 0; particleIndexRow < rowcolSize; ++particleIndexRow)
  {
	  float stepLength = 1.0f/rowcolSize;
	  for(int particleIndexCol = 0; particleIndexCol < rowcolSize; ++particleIndexCol)
	  {
		  particles[particleIndexCol + rowcolSize*particleIndexRow].m_x = -kOffset + particleIndexCol*stepLength; // + pos_dist(eng);
		  particles[particleIndexCol + rowcolSize*particleIndexRow].m_y = -kOffset + particleIndexRow*stepLength; // + pos_dist(eng);
	  }
  }
}

void borderParticlesInit()
{
	float stepLengthx = 2.0f/kBorderSideParticlesCount;
	float stepLengthy = 2.0f/kBorderSideParticlesCount;

	float mass = 0.1f;
	float md = 100.0f;
	float pressure = 0.01f;

	for(int i = 0; i < kBorderSideParticlesCount; i++)
	{
		borderParticles[i].m_x = -1.0f + stepLengthx*i;
		borderParticles[i].m_y = -0.98f;
		borderParticles[i].m_mass = mass;
		borderParticles[i].m_massDensity = md;
		borderParticles[i].m_pressure = pressure;
	}

	for(int i = 0; i < kBorderSideParticlesCount; i++)
	{
		borderParticles[kBorderSideParticlesCount + i].m_x = -0.98f;
		borderParticles[kBorderSideParticlesCount + i].m_y = (-1.0f) + stepLengthy*i;
		borderParticles[kBorderSideParticlesCount + i].m_mass = mass;
		borderParticles[kBorderSideParticlesCount + i].m_massDensity = md;
		borderParticles[kBorderSideParticlesCount + i].m_pressure = pressure;
	}

	for(int i = 0; i < kBorderSideParticlesCount; i++)
	{
		borderParticles[2*kBorderSideParticlesCount + i].m_x = 0.98f;
		borderParticles[2*kBorderSideParticlesCount + i].m_y = -1.0f + stepLengthy*i;
		borderParticles[2*kBorderSideParticlesCount + i].m_mass = mass;
		borderParticles[2*kBorderSideParticlesCount + i].m_massDensity = md;
		borderParticles[2*kBorderSideParticlesCount + i].m_pressure = pressure;
	}

	for(int i = 0; i < kBorderSideParticlesCount; i++)
	{
		borderParticles[3*kBorderSideParticlesCount + i].m_x = -1.0f + stepLengthx*i;
		borderParticles[3*kBorderSideParticlesCount + i].m_y = 0.98f;
		borderParticles[3*kBorderSideParticlesCount + i].m_mass = mass;
		borderParticles[3*kBorderSideParticlesCount + i].m_massDensity = md;
		borderParticles[3*kBorderSideParticlesCount + i].m_pressure = pressure;
	}

}

void updateGrid()
{
	memset(&grid[0], 0, kGridCellCount*sizeof(particle*));
	//grid.swap( std::vector<particle*>(grid.size(), 0) );
	for(size_t i = 0; i < kParticlesCount; i++)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;


		if (x < 1)
			x = 1;
		else if (x > kGridWidth-2)
			x = kGridWidth-2;

		if (y < 1)
			y = 1;
		else if (y > kGridHeight-2)
			y = kGridHeight-2;

		pi.next = grid[x+y*kGridWidth];
		grid[x+y*kGridWidth] = &pi;

		gridCoords[i*2] = x;
		gridCoords[i*2+1] = y;
	}
}

void updatNeighbours()
{
	for(int i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*kGridWidth;
		
		neighbours[i].count = 0;
		
		//Loop over border
		for(size_t j = 0; j < kBorderParticlesCount; j++)
		{
			particle& bp = borderParticles[j];
			float pm = bp.m_mass;

			float dx = pi.m_x - bp.m_x;
			float dy = pi.m_y - bp.m_y;
			float distance2 = dx*dx + dy*dy;

			if(distance2 < IR2)
			{
				if(neighbours[i].count < kMaxNeighbourCount)
				{
					neighbours[i].particles[neighbours[i].count] = &bp;
					neighbours[i].r2[neighbours[i].count] = distance2;
					++neighbours[i].count;
					//std::cout << "I'm on the border" << std::endl;
				}
			}
		}
		//loop over cells
		for (int ni=gi-1; ni<=gi+1; ++ni)
		{
			for (int nj=gj-kGridWidth; nj<=gj+kGridWidth; nj+=kGridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					//do fancy math
					//std::cout << "ppj x: " << ppj->m_x << std::endl;
					dx = pi.m_x - ppj->m_x;
					dy = pi.m_y - ppj->m_y;
					distance2 = dx*dx + dy*dy;

					if(distance2 < IR2)
					{
						//Density
						//massDensity += particleMass*kWdeafult* (IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2);

						if(neighbours[i].count < kMaxNeighbourCount)
						{
							neighbours[i].particles[neighbours[i].count] = ppj;
							neighbours[i].r2[neighbours[i].count] =distance2;
							++neighbours[i].count;
						}
					}
				}
			}
		}
	}
}


float calculateMass()
{
	float density = 0.0f; 
	

	for(int i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*kGridWidth;

		//loop over cells 
		for (size_t ni=gi-1; ni<=gi+1; ++ni)
		{
			for (size_t nj=gj-kGridWidth; nj<=gj+kGridWidth; nj+=kGridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					const particle& pj = *ppj;
					dx = pi.m_x - pj.m_x;
					dy = pi.m_y - pj.m_y;
					distance2 = dx*dx + dy*dy;
					if(distance2 < IR2)
					{
						density += (IR2 - distance2)*(IR2 - distance2)*(IR2 - distance2)*kWdeafult;
					}
				}
			}
		}
	}
	float dA = density/kParticlesCount;
	float mass = (dA*restDensity)/(dA*dA);
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];
		pi.m_mass = mass;
	}
	return mass;
}
