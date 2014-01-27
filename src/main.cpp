
#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <glfw3.h> // GLFW helper library
#include <stdio.h>
#include <iostream>
#include <random>

#include "particles.h"
#include "shader.h"

#define kParticlesCount 1000
#define kWindowWidth 800
#define kWindowHeight 600
float kDt = 0.001f;

void advance();
void render();
void glInit();
void particlesInit();

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

particleSystem pSys;

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
    pSys[i].update(kDt);

  }
    //std::cout << "x: " << pSys[0].m_x << " y: " << pSys[0].m_y << " u: " << pSys[0].m_u << "  v: " << pSys[0].m_v << std::endl;

}

void render()
{
    
	vbo = 0;
    glGenBuffers (1, &vbo);
    glBindBuffer (GL_ARRAY_BUFFER, vbo);
    glBufferData (GL_ARRAY_BUFFER, kParticlesCount * sizeof(particle), &pSys[0], GL_STATIC_DRAW);

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

    glPointSize(1.0f);
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

  //give all particles random starting positions and velocities
  for(int i = 0; i < kParticlesCount; i++)
  {
    float x = pos_dist(eng);
    float y = pos_dist(eng);
    float u = vel_dist(eng);
    float v = vel_dist(eng);

    particle p = particle(x,y,u,v);
    pSys.push_back(p);
    //std::cout << "x: " << p.m_x << " y: " << p.m_y << " u: " << p.m_u << "  v: " << p.m_v << std::endl;
  }
  
}