
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>

#include "particles.h"
#include "shader.h"

#include <vector>
#include <math.h>

#define kParticlesCount 1024
#define kWindowWidth 640
#define kWindowHeight 480
#define kPi 3.14159265359
#define g -9.81


void advance();
void render();
void glInit();
void particlesInit();
void drawGrid();
void updateGrid();
void createDrawablePoints();
void loopStructure();

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
const float interactionRadius =  0.1f;
const float cellSize = interactionRadius;
const float kDt = 0.001f;
const int kCellCount = 100;
const float restDensity = 988.0f;
const int kstiffnes = 20;
const float surfaceTension = 0.0728f;
const float viscosityConstant = 3.5f;
const float damp = 0.2f;

float particleMass;

float surfaceLimit = 7.0f;	// defined as sqrt(restDensity/x) 
							// where x is average particle sum in kernel
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

	particle* next;
};


struct point
{
	float x;
	float y;
};

typedef std::vector<particle> pVec;
pVec pSys;

particle particles[kParticlesCount];
//particle* grid[kCellCount][kCellCount];

const size_t kGridWidth = (size_t)(2.0 / cellSize);
const size_t kGridHeight = (size_t)(2.0 / cellSize);
const size_t kGridCellCount = kGridWidth * kGridHeight;
size_t gridCoords[kParticlesCount*2];
particle* grid[400];



std::vector<point> drawablePoints;

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

  particlesInit();

  updateGrid();
  particleMass = calculateMass();
  std::cout << "Mass: "<< particleMass << std::endl;

  glInit();


  while (!glfwWindowShouldClose (window)) 
  {
	  updateGrid();
	  loopStructure();

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
    GLuint programID = LoadShader( "default.vert", "flat.frag" );
    glUseProgram (programID);
}



void createDrawablePoints()
{
	drawablePoints.clear();
	drawablePoints.reserve(kParticlesCount);
	for(int i = 0; i < kParticlesCount; ++i)
	{
		point p;
		p.x = particles[i].m_x;
		p.y = particles[i].m_y;
		drawablePoints[i] = p;
	}
}

void render()
{
    
	vbo = 0;
    glGenBuffers (1, &vbo);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glBufferData (GL_ARRAY_BUFFER, kParticlesCount * sizeof(point), &drawablePoints[0], GL_STATIC_DRAW);

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
    glDrawArrays (GL_POINTS, 0, kParticlesCount);
}


void particlesInit()
{

  for(int particleIndexRow = 0; particleIndexRow < 32; ++particleIndexRow)
  {
	  float stepLength = 1.0f/32.0f;
	  for(int particleIndexCol = 0; particleIndexCol < 32; ++particleIndexCol)
	  {
		  particles[particleIndexCol + 32*particleIndexRow].m_x = -0.98f + particleIndexCol*stepLength;
		  particles[particleIndexCol + 32*particleIndexRow].m_y = -0.98f + particleIndexRow*stepLength;
	  }
  }
}

void updateGrid()
{
	memset(grid, 0, 400*sizeof(particle*));

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


void loopStructure()
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
					}
				}
			}
		}
		//save massDensity
		pi.m_massDensity = massDensity;
		//save Pressure
		pi.m_pressure = kstiffnes * (massDensity - restDensity);
	}


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
		//loop over cells 
		for (size_t ni=gi-1; ni<=gi+1; ++ni)
		{
			for (size_t nj=gj-kGridWidth; nj<=gj+kGridWidth; nj+=kGridWidth)
			{
				//loop over neighbors
				for (particle* ppj=grid[ni+nj]; NULL!=ppj; ppj=ppj->next)
				{
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < interactionRadius*interactionRadius)
					{
						float* Wnormal = WgradDefult(dx, dy);
						mdj = ppj->m_massDensity;
						//std::cout <<"W[0] = " << Wnormal[0] << std::endl;
						//std::cout <<"W[1] = " << Wnormal[1] << std::endl;

						normalx += (particleMass/mdj)*Wnormal[0];
						normaly += (particleMass/mdj)*Wnormal[1];

						gradNormal += (particleMass/mdj)*WlaplacianDefult(distance2);
						
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
							
							viscosityForcex += velocityDiffu * (particleMass/mdj) * WlaplacianViscosity(distance2);
							viscosityForcey += velocityDiffv * (particleMass/mdj) * WlaplacianViscosity(distance2);
							
						}
					}
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
		//Time integration
		//Sebastian Superior integration
		
		/*if(t == 0)
		{
			pi.m_u += pi.m_u + 0.5*kDt*accelerationX;
			pi.m_v += pi.m_v + 0.5*kDt*accelerationY;

			pi.m_x += kDt*pi.m_u;
			pi.m_y += kDt*pi.m_v;
		}else
		{*/
			pi.m_u += kDt*accelerationX;
			pi.m_v += kDt*accelerationY;

			pi.m_x += kDt*pi.m_u;
			pi.m_y += kDt*pi.m_v;
		//}

		//Colision handling and response
		float current,cp,d,n,u;
		if(pi.m_x < -1)
		{
			current = pi.m_x;
			u = pi.m_u;
			cp = -1;

			d = sqrt((cp-current)*(cp-current));
			n = 1;
			pi.m_x = cp + d*n;
			pi.m_u = u - (1 + damp*(d/(kDt*sqrt(u*u))))*(u*n)*n;
		}

		if(pi.m_x > 1)
		{
			current = pi.m_x;
			u = pi.m_u;
			cp = 1;

			d = sqrt((cp-current)*(cp-current));
			n = -1;

			pi.m_x = cp + d*n;
			pi.m_u = u - (1 + damp*(d/(kDt*sqrt(u*u))))*(u*n)*n;
		}

		if(pi.m_y < -1)
		{
			current = pi.m_y;
			u = pi.m_v;
			cp = -1;

			d = sqrt((cp-current)*(cp-current));
			n = 1;

			pi.m_y = cp + d*n;
			pi.m_v = u - (1 + damp*(d/(kDt*sqrt(u*u))))*(u*n)*n;

		}
		
		if(pi.m_y > 1)
		{
			current = pi.m_y;
			u = pi.m_v;
			cp = 1;

			d = sqrt((cp-current)*(cp-current));
			n = -1;

			pi.m_y = cp + d*n;
			pi.m_v = u - (1 + damp*(d/(kDt*sqrt(u*u))))*(u*n)*n;
		}

	}
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
	return (dA*restDensity)/(dA*dA);
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

/*
float distance2(particle* currentParticle, particle* neighbour)
{
	// |ri - rj|^2
	return (currentParticle->m_x - neighbour->m_x)*(currentParticle->m_x - neighbour->m_x) +
		(currentParticle->m_y - neighbour->m_y)*(currentParticle->m_y - neighbour->m_y);
}
*/