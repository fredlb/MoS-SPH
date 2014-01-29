
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>

#include "particles.h"
#include "shader.h"

#include <vector>
#include <math.h>

#define kParticlesCount 1000
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
const float kParticleMass = 0.02f;
const float restDensity = 988.0f;
const int kstiffnes = 5;
const float surfaceTension = 0.07f;
const float viscosityConstant = 3.5f;

float surfaceLimit = 0.0f;	// defined as sqrt(restDensity/x) 
							// where x is average particle sum in kernel
float accelerationX;
float accelerationY;

//Bör vara en del av particle structen så att man kommer åt "partikel j"
float massDensityArray[kParticlesCount]; 
float pressureArray[kParticlesCount];


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
particle* grid[kCellCount][kCellCount];


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
  //particlesInit();
  glInit();

  while (!glfwWindowShouldClose (window)) 
  {

    advance();

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

void advance()
{
  for(int i = 0; i < kParticlesCount; i++)
  {
    //pSys[i].update(kDt);

  }
    //std::cout << "x: " << pSys[0].m_x << " y: " << pSys[0].m_y << " u: " << pSys[0].m_u << "  v: " << pSys[0].m_v << std::endl;

}

void render()
{
    
	vbo = 0;
    glGenBuffers (1, &vbo);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glBufferData (GL_ARRAY_BUFFER, kParticlesCount * sizeof(particle), &particles[0], GL_STATIC_DRAW);

    vao = 0;
    glGenVertexArrays (1, &vao);
    glBindVertexArray (vao);
    glEnableVertexAttribArray (0);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer (0, 4, GL_FLOAT, GL_FALSE, 0, NULL);

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
    glDrawArrays (GL_POINTS, 0, pSys.size());
}

/**
* just distribute particles randomly in [-1,1] with random velocities
*/
void particlesInit()
{
	//PRNG (c++11)
  std::mt19937 eng((std::random_device())());
  std::uniform_real_distribution<> pos_dist(-1,1);
  std::uniform_real_distribution<> vel_dist(-0.3,0.3);
}

void updateGrid()
{
	memset(grid, 0, kCellCount*sizeof(particle*));

	for(size_t i = 0; i < kParticlesCount; i++)
	{
		particle& pi = particles[i];

		int x = pi.m_x/cellSize;
		int y = pi.m_y/cellSize;

		if(x < 1)
		{
			x = 1;
		}
		else if( x > kCellCount - 2)
		{
			x = kCellCount - 2;
		}

		if(y < 1)
		{
			y = 1;
		}
		else if( y > kCellCount - 2)
		{
			y = kCellCount - 2;
		}

		pi.next = grid[x][y];
		grid[x][y] = &pi;
	}
}


void loopStructure()
{
	//Mass-density and pressure loop
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];
		int x = pi.m_x/cellSize;
		int y = pi.m_y/cellSize;

		float massDensity = 0.0;
		//loop over cells 
		for(size_t gx = x-1; gx < x+1; gx++)
		{
			for(size_t gy = y-1; gy < y+1; gy++)
			{
				//loop over neighbors
				for (particle* ppj=grid[gx][gy]; NULL!=ppj; ppj=ppj->next)
				{
					//do fancy math
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;

					if(distance2 < interactionRadius*interactionRadius)
					{
						//Density
						massDensity += kParticleMass*Wdeafult(distance2);
					}
				}
			}
		}
		//save massDensity
		massDensityArray[i] = massDensity;
		//save Pressure
		pressureArray[i] = kstiffnes * (massDensity - restDensity);
	}


	//Force loop
	for(size_t i = 0; i < kParticlesCount; ++i)
	{
		particle& pi = particles[i];
		int x = pi.m_x/cellSize;
		int y = pi.m_y/cellSize;

		pressureForcex = 0.0f; 
		pressureForcey = 0.0f;
		viscosityForcex = 0.0f;
		viscosityForcey = 0.0f;
		normalx = 0.0f;
		normaly = 0.0f;
		gradNormal = 0.0f;
		surfaceTensionForcex = 0.0f;
		surfaceTensionForcey = 0.0f;
		gravity = g*massDensityArray[i];

		//loop over cells 
		for(size_t gx = x-1; gx < x+1; gx++)
		{
			for(size_t gy = y-1; gy < y+1; gy++)
			{
				//loop over neighbors
				for (particle* ppj=grid[gx][gy]; NULL!=ppj; ppj=ppj->next)
				{
					float dx = pi.m_x - ppj->m_x;
					float dy = pi.m_y - ppj->m_y;
					float distance2 = dx*dx + dy*dy;
					if(distance2 < interactionRadius*interactionRadius)
					{
						float* Wnormal = WgradDefult(dx, dy);
						/*Need acces to ppj massDensity
						
						normalx += (kParticleMass/massDensityArray[j])*Wnormal[0];
						normaly += (kParticleMass/massDensityArray[j])*Wnormal[1];

						gradNormalx += (kParticleMass/massDensityArray[j]*WlaplacianDefult(distance2);
						*/
						if( distance2 != 0)
						{
							float* W = WgradPressure(dx,dy);
							float velocityDiffu = pi.m_u - ppj->m_u;
							float velocityDiffv = pi.m_v - ppj->m_v;
							/* Need acces to ppj massDensity and Pressure
							pressureForcex += ((pressureArray[i]/pow(massDensityArray[i],2))+(pressureArray[j]/pow(massDensityArray[j],2))*kparticleMass*W[0];
							pressureForcey += ((pressureArray[i]/pow(massDensityArray[i],2))+(pressureArray[j]/pow(massDensityArray[j],2))*kparticleMass*W[1];
							
							viscosityForcex += velocityDiffu * (kParticleMass/massDensityArray[j]) * WlaplacianViscosity(distance2);
							viscosityForcey += velocityDiffu * (kParticleMass/massDensityArray[j]) * WlaplacianViscosity(distance2);
							*/
						}
					}
				}
			}
		}

		float normalLenght = sqrt(normalx*normalx + normaly*normaly);
		if(normalLenght > surfaceLimit){
			surfaceTensionForcex = - surfaceTension  * gradNormal * (normalx/normalLenght);
		}
		pressureForcex = -massDensityArray[i] * pressureForcex;
		pressureForcey = -massDensityArray[i] * pressureForcey;
		
		viscosityForcex = viscosityConstant*viscosityForcex;
		viscosityForcey = viscosityConstant*viscosityForcey;
		
		accelerationX = (pressureForcex + viscosityForcex + surfaceTensionForcex)/massDensityArray[i];
		accelerationY = (pressureForcex + viscosityForcey + surfaceTensionForcey + gravity)/massDensityArray[i];
		
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
			pi.m_v += kDt*pi.m_v;
		//}

	}
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
	W[1] = -(45/(kPi*pow(interactionRadius,6)))*(dx/sqrt(distance2))*pow((interactionRadius-sqrt(distance2)),2);

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
	W[1] = -(945/(32*kPi*pow(interactionRadius,9)))*dx*pow((pow(interactionRadius,2)-distance2),2);
	
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