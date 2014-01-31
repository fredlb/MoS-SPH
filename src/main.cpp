
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>
#include <algorithm>

#include "particles.h"
#include "shader.h"

#include <vector>
#include <math.h>

#define kParticlesCount 1024
#define kBorderParticlesCount 300
#define kWindowWidth 640
#define kWindowHeight 480
#define kPi 3.14159265359
#define g -9.81
#define kFrameRate 60
#define kSubSteps 7

#define kOffset 0.5f

//#define kDt ((1.0f/kFrameRate) / kSubSteps)

#define averageParticles 20
//#define interactionRadius sqrt(averageParticles/(kParticlesCount*kPi))
#define interactionRadius 0.1f
#define cellSize (2.0f*interactionRadius)

void advance();
void render();
void glInit();
void particlesInit();
void borderParticlesInit();
void drawGrid();
void updateGrid();
void createDrawablePoints();
void calculatePressure();
void integrate();
void calulateForces();

float calculateMass();

float Wdeafult(float distance2);
float* WgradPressure(float dx, float dy);
float WlaplacianViscosity(float distance2);
float* WgradDefult(float dx, float dy);
float WlaplacianDefult(float distance2);

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

const float kViewScale =  2.0f;
//const float interactionRadius =  0.05f;
//const float cellSize = 2*interactionRadius;



const float kDt = 0.0001f;
const int kCellCount = 100;
const float restDensity = 988.0f;
const int kstiffnes = 50;
const float surfaceTension = 0.0728f;
const float viscosityConstant = 3.5f;
const float damp = 0.2f;

float particleMass;

float surfaceLimit = sqrt(restDensity/averageParticles);
float accelerationX;
float accelerationY;

//Every force
float	pressureForcex,	pressureForcey,	viscosityForcex,
		viscosityForcey,normalx, normaly, gradNormal,
		surfaceTensionForcex,surfaceTensionForcey, gravity;


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
    const particle* particles[kMaxNeighbourCount];
    float r[kMaxNeighbourCount];
    size_t count;
};

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
std::vector<particle*> grid;



std::vector<point> drawablePoints;
std::vector<point> DEBUG_CORNER;

point acceleration[kParticlesCount];
point prevAcceleration[kParticlesCount];
float vhx[kParticlesCount];
float vhy[kParticlesCount];
bool firstIteration = true;

Neighbours neighbours[kParticlesCount];


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
  
  grid.reserve(kGridCellCount);
  particlesInit();
  borderParticlesInit();
  updateGrid();

  particleMass = calculateMass();
  std::cout << "Mass: "<< particleMass << std::endl;
  

  glInit();

  std::cout << "interaction radius: " << interactionRadius << std::endl;
  std::cout << "cellSize: " << cellSize << std::endl;
  std::cout << "grid cell count: " << kGridCellCount << std::endl;
  std::cout << "kgridwidth: " << kGridWidth << std::endl;
  std::cout << "kDt: " << kDt << std::endl;


  double t = 0.0;
  double currentTime = glfwGetTime();
  double accumulator = 0.0;

  while (!glfwWindowShouldClose (window)) 
  {
	  /*double newTime = glfwGetTime();

	  double frameTime = newTime - currentTime;
	  currentTime = newTime;

	  accumulator += frameTime;*/

	  
	  //while(accumulator >= kDt)
	  for(int i = 0; i < kSubSteps; ++i)
	  {
		  updateGrid();
		  calculatePressure();
		  calulateForces();
		  //integrate();

		  //accumulator -= kDt;
		  //t += kDt;
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
    GLuint programID = LoadShader( "src/default.vert", "src/flat.frag" );
    glUseProgram (programID);


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




}



void createDrawablePoints()
{
	drawablePoints.clear();
	drawablePoints.reserve(kParticlesCount + kBorderParticlesCount);
	for(int i = 0; i < kParticlesCount; ++i)
	{
		point p;
		p.x = particles[i].m_x;
		p.y = particles[i].m_y;
		drawablePoints[i] = p;
	}
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

    
	vbo = 0;
    glGenBuffers (1, &vbo);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glBufferData (GL_ARRAY_BUFFER, (kParticlesCount + kBorderParticlesCount) * sizeof(point), &drawablePoints[0], GL_STATIC_DRAW);

    vao = 0;
    glGenVertexArrays (1, &vao);
    glBindVertexArray (vao);
    glEnableVertexAttribArray (0);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer (0, 2, GL_FLOAT, GL_FALSE, 0, NULL);

    glClear (GL_COLOR_BUFFER_BIT);
    glBindVertexArray (vao);

    glClearColor(0.05f, 0.05f, 0.05f, 1);
 
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
}


void particlesInit()
{

	std::mt19937 eng((std::random_device())());
	std::uniform_real_distribution<> pos_dist(-0.001,0.001);

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
	float stepLengthx = 2.0f/kBorderParticlesCount;
	float stepLengthy = interactionRadius/4;
	for(int i = 0; i<kBorderParticlesCount; i++)
	{
	//for(int j = 0; j<4; j++)
	//{
		borderParticles[i].m_x = -1.0f + stepLengthx*i;
		borderParticles[i].m_y = -0.98f;
		borderParticles[i].m_mass = 1;
		borderParticles[i].m_massDensity = 10*restDensity;
		borderParticles[i].m_pressure = kstiffnes*(9*restDensity);

	//}
	}
}

void updateGrid()
{
	memset(&grid[0], 0, kGridCellCount*sizeof(particle*));

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

void calulateForces()
{
	//Force loop
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		pressureForcex = 0.0f; 
		pressureForcey = 0.0f;
		viscosityForcex = 0.0f;
		viscosityForcey = 0.0f;
		normalx = 0.0f;
		normaly = 0.0f;
		gradNormal = 0.0f;
		surfaceTensionForcex = 0.0f;
		surfaceTensionForcey = 0.0f;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*kGridWidth;

		float mdi = pi.m_massDensity;
		float mdj = 0.0f;
		gravity = g*mdi;
		//loop over neighbours 
				for(size_t j=0; j < neighbours[i].count; ++j)
				{
					const particle* ppj = neighbours[i].particles[j];
					float massi = pi.m_mass;
					float massj = ppj->m_mass;
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < interactionRadius*interactionRadius)
					{
						float* Wnormal = WgradDefult(dx, dy);
						mdj = ppj->m_massDensity;

						normalx += (massj/mdj)*Wnormal[0];
						normaly += (massj/mdj)*Wnormal[1];

						gradNormal += (massj/mdj)*WlaplacianDefult(distance2);

						if( distance2 != 0)
						{
							float* W = WgradPressure(dx,dy);
							// uj - ui
							float velocityDiffu = ppj->m_u - pi.m_u;
							float velocityDiffv = ppj->m_v - pi.m_v;

							float pressi = pi.m_pressure;
							float pressj = ppj->m_pressure;
							pressureForcex += ((pressi/pow(mdi,2))+(pressj/pow(mdj,2)))*particleMass*W[0];
							pressureForcey += ((pressi/pow(mdi,2))+(pressj/pow(mdj,2)))*particleMass*W[1];

							viscosityForcex += velocityDiffu * (massj/mdj) * WlaplacianViscosity(distance2);
							viscosityForcey += velocityDiffv * (massj/mdj) * WlaplacianViscosity(distance2);

						}
					}
				}
		float normalLenght = sqrt(normalx*normalx + normaly*normaly);
		if(normalLenght > surfaceLimit){
			surfaceTensionForcex = - surfaceTension  * gradNormal * (normalx/normalLenght);
			surfaceTensionForcey = - surfaceTension  * gradNormal * (normaly/normalLenght);
		}
		pressureForcex = -mdi * pressureForcex;
		pressureForcey = -mdi * pressureForcey;

		viscosityForcex = viscosityConstant*viscosityForcex;
		viscosityForcey = viscosityConstant*viscosityForcey;

		accelerationX = (pressureForcex + viscosityForcex + surfaceTensionForcex)/mdi;
		accelerationY = (pressureForcey + viscosityForcey + surfaceTensionForcey + gravity)/mdi;

		//True leap-frog
                if(firstIteration){
                        vhx[i] = pi.m_u + 0.5*accelerationX*kDt;
                        vhy[i] = pi.m_v + 0.5*accelerationY*kDt;
                        
                        pi.m_u += accelerationX*kDt;
                        pi.m_v += accelerationY*kDt;

                        pi.m_x += vhx[i]*kDt;
                        pi.m_y += vhy[i]*kDt;

                        firstIteration = false;
                }else{
                        vhx[i] += accelerationX*kDt;
                        vhy[i] += accelerationY*kDt;
                        
                        pi.m_u = vhx[i] + 0.5*accelerationX*kDt;
                        pi.m_v = vhy[i] + 0.5*accelerationY*kDt;

                        pi.m_x += vhx[i]*kDt;
                        pi.m_y += vhy[i]*kDt;
                }

	}

		
}

void calculatePressure()
{
	//Mass-density and pressure loop
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*kGridWidth;
		
		float massDensity = 0.0f;
		neighbours[i].count = 0;
		
		//Loop over border
		for(size_t j = 0; j < kBorderParticlesCount; j++)
		{
			particle& bp = borderParticles[j];
			float pm = bp.m_mass;

			float dx = pi.m_x - bp.m_x;
			float dy = pi.m_y - bp.m_y;
			float distance2 = dx*dx + dy*dy;

			if(distance2 < interactionRadius*interactionRadius)
			{
				//Density
				massDensity += particleMass*Wdeafult(distance2);

				if(neighbours[i].count < kMaxNeighbourCount)
				{
					neighbours[i].particles[neighbours[i].count] = &bp;
					neighbours[i].r[neighbours[i].count] = sqrt(distance2);
					++neighbours[i].count;
					//std::cout << "I'm on the border" << std::endl;
				}
			}
		}
		//loop over cells 
		for (size_t ni=gi-1; ni<=gi+1; ++ni)
		{
			for (size_t nj=gj-kGridWidth; nj<=gj+kGridWidth; nj+=kGridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					//do fancy math
					//std::cout << "ppj x: " << ppj->m_x << std::endl;
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;

					if(distance2 < interactionRadius*interactionRadius)
					{
						//Density
						massDensity += particleMass*Wdeafult(distance2);

						if(neighbours[i].count < kMaxNeighbourCount)
						{
							neighbours[i].particles[neighbours[i].count] = ppj;
							neighbours[i].r[neighbours[i].count] = sqrt(distance2);
							++neighbours[i].count;
						}
					}
				}
			}
		}
		//save massDensity
		pi.m_massDensity = massDensity;
		//std::cout << massDensity << std::endl;
		//save Pressure
		pi.m_pressure = kstiffnes * (massDensity - restDensity);
		//std::cout << pi.m_pressure << std::endl;
	}


	
}

void integrate()
{
	for (int i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		//True leap-frog
		if(firstIteration){
			vhx[i] = pi.m_u + 0.5*accelerationX*kDt;
			vhy[i] = pi.m_v + 0.5*accelerationY*kDt;

			pi.m_u += accelerationX*kDt;
			pi.m_v += accelerationY*kDt;

			pi.m_x += vhx[i]*kDt;
			pi.m_y += vhy[i]*kDt;

			firstIteration = false;
		}else{
			vhx[i] += accelerationX*kDt;
			vhy[i] += accelerationY*kDt;

			pi.m_u = vhx[i] + 0.5*accelerationX*kDt;
			pi.m_v = vhy[i] + 0.5*accelerationY*kDt;

			pi.m_x += vhx[i]*kDt;
			pi.m_y += vhy[i]*kDt;
		}

		/*
		pi.m_x += pi.m_u*kDt + 0.5*acceleration[i].x*kDt*kDt;
		pi.m_y += pi.m_v*kDt + 0.5*acceleration[i].y*kDt*kDt;


		pi.m_u += 0.5*(acceleration[i].x + prevAcceleration[i].x)*kDt;
		pi.m_v += 0.5*(acceleration[i].y + prevAcceleration[i].y)*kDt;

		prevAcceleration[i].x = acceleration[i].x;
		prevAcceleration[i].y = acceleration[i].y;
		*/
		//Colision handling and response
		float current,cp,d,n,u,v;
		if(pi.m_x < -1)
		{
			current = pi.m_x;
			u = pi.m_u;
			v = pi.m_v;
			cp = -1;

			d = sqrt((cp-current)*(cp-current));
			n = 1;
			pi.m_x = cp + d*n;
			vhx[i] = u - (1 + damp*(d/(kDt*sqrt(u*u+v*v))))*(u*n)*n;
			vhy[i] = v;
		}

		if(pi.m_x > 1)
		{
			current = pi.m_x;
			u = pi.m_u;
			v = pi.m_v;
			cp = 1;

			d = sqrt((cp-current)*(cp-current));
			n = -1;

			pi.m_x = cp + d*n;
			vhx[i] = u - (1 + damp*(d/(kDt*sqrt(u*u+v*v))))*(u*n)*n;
			vhy[i] = v;
		}
		
		if(pi.m_y < -1)
		{
			current = pi.m_y;
			v = pi.m_v;
			u = pi.m_u;
			cp = -1;

			d = sqrt((cp-current)*(cp-current));
			n = 1;

			pi.m_y = cp + d*n;
			vhy[i] = v - (1 + damp*(d/(kDt*sqrt(u*u+v*v))))*(v*n)*n;
			vhx[i] = u;

		}
		
		if(pi.m_y > 1)
		{
			current = pi.m_y;
			u = pi.m_u;
			v = pi.m_v;
			cp = 1;

			d = sqrt((cp-current)*(cp-current));
			n = -1;

			pi.m_y = cp + d*n;
			vhy[i] = v - (1 + damp*(d/(kDt*sqrt(u*u+v*v))))*(v*n)*n;
			vhx[i] = u;
		}
		
	}
		
}

void collisionResponse()
{

}

float calculateMass()
{
	float density = 0.0f; 
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];

		int x = (1 + pi.m_x)/cellSize;
		int y = (1 + pi.m_y)/cellSize;

		size_t gi = gridCoords[i*2];
		size_t gj = gridCoords[i*2+1]*kGridWidth;
		
		/* For Fredriks inferior integration
		prevAcceleration[i].x = 0.0f;
		prevAcceleration[i].y = 0.0f;

		acceleration[i].x = 0.0f;
		acceleration[i].y = 0.0f;
		*/

		//loop over cells 
		for (size_t ni=gi-1; ni<=gi+1; ++ni)
		{
			for (size_t nj=gj-kGridWidth; nj<=gj+kGridWidth; nj+=kGridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					const particle& pj = *ppj;
					//std::cout << "pi x:" << pi.m_x << std::endl;

					//std::cout << "pj x:" << pj.m_x << std::endl;
					float dx = pi.m_x - pj.m_x;
					//std::cout << "dx: " << dx << std::endl;
					float dy = pi.m_y - pj.m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < interactionRadius*interactionRadius)
					{
						density += Wdeafult(distance2);
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

float WKernel(float distance2)
{
	return 0.0f;
}

float Wdeafult(float distance2)
{

	float W = (315/(64*kPi*pow(interactionRadius,9))) * pow((pow(interactionRadius,2) -distance2),3);
	return W;
}

float* WgradPressure(float dx, float dy)
{
	float W[2];
	float distance2 = dx*dx + dy*dy;
	
	W[0] = -(45/(kPi*pow(interactionRadius,6)))*(dx/sqrt(distance2))*pow((interactionRadius-sqrt(distance2)),2);
	W[1] = -(45/(kPi*pow(interactionRadius,6)))*(dy/sqrt(distance2))*pow((interactionRadius-sqrt(distance2)),2);

	return W;
}

float WlaplacianViscosity(float distance2)
{
	float W = (45/(kPi*pow(interactionRadius,6)))*(interactionRadius-sqrt(distance2));
	return W;
}

float* WgradDefult(float dx, float dy)
{
	float W[2];
	float distance2 = dx*dx + dy*dy;

	W[0] = -(945/(32*kPi*pow(interactionRadius,9)))*dx*pow((pow(interactionRadius,2)-distance2),2);
	W[1] = -(945/(32*kPi*pow(interactionRadius,9)))*dy*pow((pow(interactionRadius,2)-distance2),2);
	
	return W;
}

float WlaplacianDefult(float distance2)
{
	float W = -(945/(32*kPi*pow(interactionRadius,9)))*(pow(interactionRadius,2)-distance2)*(3*pow(interactionRadius,2)-7*distance2);
	return W;
}
