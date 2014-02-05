#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GL/freeglut.h>

#include "shader.h"
#include "ParticleSystem.h"
#include "vec2.h"
#include "Postprocess.h"

#define kWindowWidth 640
#define kWindowHeight 480

void drawPoints(std::vector<vec2> points, float r, float g, float b, float a, float size);
void drawMetaballs(std::vector<vec2> points);
void render();
void update();
void glInit();

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

ParticleSystem* simulation;
Postprocess* metaballprocess;

const float kViewScale =  2.0f;

GLuint programID = 0;
GLuint gradientProgramID = 0;
std::vector<vec2> pointsToDraw;
int current_time;

int main (int argc, char** argv) {
	current_time = glutGet(GLUT_ELAPSED_TIME);
	simulation = new ParticleSystem();
	//for(int i=0; i<1; i++) simulation->advance();
	glutInitWindowSize(kWindowWidth, kWindowHeight);
	glutInit(&argc, argv);
	glutInitDisplayString("samples stencil>=3 rgb double depth");
	glutCreateWindow("SPH");
	glutDisplayFunc(render);
	glutIdleFunc(update);
               
	// start GLEW extension handler
	glewExperimental = GL_TRUE;
	glewInit ();

	// get version info
	const GLubyte* renderer = glGetString (GL_RENDERER); // get renderer string
	const GLubyte* version = glGetString (GL_VERSION); // version as a string
	printf ("Renderer: %s\n", renderer);
	printf ("OpenGL version supported %s\n", version);



	glInit();

	glutMainLoop();

	return 0;
}

void glInit()
{
    programID = LoadShader( "src/default.vert", "src/flat.frag" );
	gradientProgramID = LoadShader( "src/default.vert", "src/gradientball.frag" );
	vbo = 0;
	glGenBuffers (1, &vbo);
	vao = 0;
	glGenVertexArrays (1, &vao);

	metaballprocess = new Postprocess(kWindowWidth, kWindowHeight, "src/metaballs.frag");
}


void render()
{
	std::string s = "Time per frame: " + std::to_string(glutGet(GLUT_ELAPSED_TIME) - current_time);
	glutSetWindowTitle(s.c_str());
	current_time = glutGet(GLUT_ELAPSED_TIME);
	//glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
	//glClear (GL_COLOR_BUFFER_BIT);
	//draw points here


	pointsToDraw = simulation->getParticleCoordinates();
	drawMetaballs(pointsToDraw);
	//drawPoints(pointsToDraw, 0.5, 0.2, 0.5, 1.0, 5.0);
	//drawPoints(simulation->getParticleCoordinatesPressure(PRESSURE_UNDER,150.0), 1.0, 0.0, 0.0, 0.5, 5.0);
	//drawPoints(simulation->getParticleCoordinatesPressure(PRESSURE_OVER,150.0), 0.0, 1.0, 0.0, 0.5, 5.0);

	//stop drawing here
	glutSwapBuffers();
}

void update()
{
	for(int i = 0; i < 3; ++i)
	{
		simulation->advance();
	}
	glutPostRedisplay();
}

void drawPoints(std::vector<vec2> points, float r, float g, float b, float a, float size)
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glUseProgram (programID);
	GLint loc = glGetUniformLocation(programID, "uColor");
    glUniform4f(loc, r, g, b, a);

	glBindBuffer (GL_ARRAY_BUFFER, vbo);
	glBufferData (GL_ARRAY_BUFFER, (points.size()) * sizeof(vec2), &points[0], GL_STATIC_DRAW);
	glBindVertexArray (vao);
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

    glPointSize(size);

    // draw points from the currently bound VAO with current in-use shader
    glDrawArrays (GL_POINTS, 0, points.size());
}

void drawMetaballs(std::vector<vec2> points)
{
	//draw to texture
	metaballprocess->drawTo();
	glClearColor(0.0, 0.0, 0.0, 1.0);
	
	glEnable(GL_POINT_SPRITE);
	glUseProgram (gradientProgramID);

	/*
	glUseProgram (programID);
	GLint loc = glGetUniformLocation(programID, "uColor");
    glUniform4f(loc, 0.5, 0.5, 1.0, 1.0);
	*/

	glBindBuffer (GL_ARRAY_BUFFER, vbo);
	glBufferData (GL_ARRAY_BUFFER, (points.size()) * sizeof(vec2), &points[0], GL_STATIC_DRAW);
	glBindVertexArray (vao);
	glEnableVertexAttribArray (0);
	glBindBuffer (GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer (0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    
 
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, 640, 0, 480, 0, 1);
 
	
	//Draw points as smooth balls (with AA)
	glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);

    glPointSize(35.0);

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    // draw points from the currently bound VAO with current in-use shader
    glDrawArrays (GL_POINTS, 0, points.size());
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//draw to screen
	metaballprocess->renderFrom();
}