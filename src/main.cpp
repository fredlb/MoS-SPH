#include <stdlib.h>
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>
#include <algorithm>
#include <vector>
#include <math.h>

#include "shader.h"
#include "ParticleSystem.h"

#define kWindowWidth 640
#define kWindowHeight 480

#define kParticlesCount 1024
#define kBorderParticlesCount 300
#define kPi 3.14159265359
#define g -9.81
#define kFrameRate 60
#define kSubSteps 7

#define kOffset 0.5f

#define averageParticles 20
//#define interactionRadius sqrt(averageParticles/(kParticlesCount*kPi))
#define interactionRadius 0.035f
#define IR2 interactionRadius*interactionRadius
#define cellSize (2.0f*interactionRadius)

void render();
void glInit();
//void drawGrid();
void createDrawablePoints(ParticleSystem s);

void advance();
void particlesInit();
void borderParticlesInit();
void updateGrid();
void calculatePressure();
void calulateForces();

float calculateMass();

float Wdeafult(float distance2);
float* WgradPressure(float dx, float dy);
float WlaplacianViscosity(float distance2);
float* WgradDefult(float dx, float dy);
float WlaplacianDefult(float distance2);

const float kWdeafult = (315/(64*kPi*pow(interactionRadius,9))); 
const float kWgradPressure = -(45/(kPi*pow(interactionRadius,6)));
const float kWlaplacianViscosity = (45/(kPi*pow(interactionRadius,6)));
const float kWgradDefult = -(945/(32*kPi*pow(interactionRadius,9)));
const float kWlaplacianDefult = -(945/(32*kPi*pow(interactionRadius,9)));

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

const float kViewScale =  2.0f;
//const float interactionRadius =  0.05f;
//const float cellSize = 2*interactionRadius;



const float kDt = 0.00025f;
const int kCellCount = 100;
const float restDensity = 988.0f;
const int kstiffnes = 100;
const float surfaceTension = 0.0728f;
const float viscosityConstant = 3.5f;
const float damp = 0.2f;

float particleMass;

float surfaceLimit = sqrt(restDensity/averageParticles);
float accelerationX;
float accelerationY;




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


bool firstIteration = true;

GLuint programID = 0;

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
  
  ParticleSystem simulation;


  drawablePoints.resize(kParticlesCount + kBorderParticlesCount);



  glInit();



  double t = 0.0;
  double currentTime = glfwGetTime();
  double accumulator = 0.0;

  while (!glfwWindowShouldClose (window)) 
  {

	  for(int i = 0; i < kSubSteps; ++i)
	  {

		  simulation.advance();

		  //integrate();

		  //accumulator -= kDt;
		  //t += kDt;
	  }

	  createDrawablePoints(simulation);
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
     programID = LoadShader( "src/default.vert", "src/flat.frag" );
    


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



void createDrawablePoints(ParticleSystem s)
{

	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < kParticlesCount; ++i)
	{
		point p;
		p.x = s.particles[i].m_x;
		p.y = s.particles[i].m_y;
		drawablePoints[i] = p;
	}
	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < kBorderParticlesCount; ++i)
	{
		point p;
		p.x = s.borderParticles[i].m_x;
		p.y = s.borderParticles[i].m_y;
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


