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

#define kWindowWidth 800
#define kWindowHeight 600

void drawPoints(std::vector<vec2> points, float r, float g, float b, float a, float size);
void drawMetaballs(std::vector<vec2> points);
void render();
void update();
void glInit();
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void keyPressed(unsigned char c, int x, int y);

unsigned int vao;
unsigned int vbo;
unsigned int shader_programme;

ParticleSystem* simulation;
Postprocess* metaballprocess;

const float kViewScale =  2.0f;

GLuint programID = 0;
GLuint gradientProgramID = 0;
std::vector<vec2> pointsToDraw;
std::vector<vec2> bordersToDraw;
int current_time;
bool POINTS_MODE = true;
int POINTS_RELOAD = true;
char keyPress = ' ';

bool lbuttonDown = false;
bool rbuttonDown = false;
int mouseX;
int mouseY;


//Controller Window
void render2win();
void valueControllerInit();
void valueControllerMouse(int button, int state, int x, int y);
void valueControllerMotion(int x, int y);
int stiffnessValueChanger = 150;
bool stiffnessValueChanged = false;
void drawGLString(GLfloat x, GLfloat y, char *textString);
int viscosityValueChanger = 150;
bool viscosityValueChanged = false;



int main (int argc, char** argv) {

	glutInit(&argc, argv);

	//controller window
	glutInitWindowSize(500, 300);
	glutCreateWindow("Value controll window");
	glutPositionWindow(1000,100);
	valueControllerInit();
	glutDisplayFunc(render2win);
	glutMouseFunc(valueControllerMouse);

	current_time = glutGet(GLUT_ELAPSED_TIME);
	simulation = new ParticleSystem();
	//for(int i=0; i<1; i++) simulation->advance();
	glutInitWindowSize(kWindowWidth, kWindowHeight);
	
	glutInitDisplayString("samples stencil>=3 rgb double depth");
	glutCreateWindow("SPH");
	glutDisplayFunc(render);
	glutIdleFunc(update);
    glutKeyboardFunc(keyPressed);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
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

void render2win(){

	glClear(GL_COLOR_BUFFER_BIT);  

	//STIFFNESS CONTROLLER
    int currentStiffness = stiffnessValueChanger*0.167;
	char buffer[5];
	sprintf(buffer,"%d", currentStiffness);
	drawGLString(70.0, 480.0, "Stiffness");
	drawGLString(154.0, 480.0, buffer);

    glBegin(GL_LINES);
	//glColor3f(1.0,1.0,0.0); 
    glVertex2d(150, 480);
    glVertex2d(150, 20);
    glEnd();

	if(stiffnessValueChanged){
		simulation->setStiffness(currentStiffness);
		stiffnessValueChanged=false;
	}

	glBegin(GL_LINES);
	glVertex2d(130, stiffnessValueChanger*1.67);
	glVertex2d(170, stiffnessValueChanger*1.67);
	glEnd();
	//END OF STIFFNESS CONTROLLER

	//VISCOSITY CONTROLLER
	int currentViscosity = viscosityValueChanger*0.033;
	//buffer.clear();
	sprintf(buffer,"%d", currentViscosity);
	drawGLString(270.0, 480.0, "Viscosity");
	drawGLString(354.0, 480.0, buffer);

    glBegin(GL_LINES);
	//glColor3f(1.0,1.0,0.0); 
    glVertex2d(350, 480);
    glVertex2d(350, 20);
    glEnd();

	if(viscosityValueChanged){
		simulation->setViscosity(currentViscosity);
		viscosityValueChanged=false;
	}

	glBegin(GL_LINES);
	glVertex2d(330, viscosityValueChanger*1.67);
	glVertex2d(370, viscosityValueChanger*1.67);
	glEnd();
	//END OF VISCOSITY CONTROLLER


	glFlush();
}

void render()
{
	std::string s = "Time per frame: " + std::to_string(glutGet(GLUT_ELAPSED_TIME) - current_time);
	glutSetWindowTitle(s.c_str());
	current_time = glutGet(GLUT_ELAPSED_TIME);
	//glClearColor(0.05f, 0.05f, 0.05f, 1.0f);
	//glClear (GL_COLOR_BUFFER_BIT);
	//draw points here


	if(!POINTS_RELOAD)
	{
		simulation->reloadParticleSystem(keyPress);
		POINTS_RELOAD = true;
	}

	pointsToDraw = simulation->getParticleCoordinates();
	bordersToDraw = simulation->getParticleCoordinatesBorder();

	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	if(POINTS_MODE)
	{
		drawPoints(pointsToDraw, 0.2, 0.55, 0.75, 1.0, 1.0);
	}
	else
	{
		drawMetaballs(pointsToDraw);
	}

	drawPoints(bordersToDraw, 0.8, 0.3, 0.3, 1.0, 2.0);

	/*
	pointsToDraw = simulation->getParticleCoordinatesBorder();
	if(POINTS_MODE)
	{
		drawPoints(pointsToDraw, 0.0, 0.55, 0.75, 1.0, 1.0);
	}
	else
	{
		drawMetaballs(pointsToDraw);
	}*/
	
	//drawPoints(simulation->getParticleCoordinatesPressure(PRESSURE_UNDER,150.0), 1.0, 0.0, 0.0, 0.5, 5.0);
	//drawPoints(simulation->getParticleCoordinatesPressure(PRESSURE_OVER,150.0), 0.0, 1.0, 0.0, 0.5, 5.0);

	//stop drawing here
	glutSwapBuffers();
}

void update()
{
	simulation->updateMouseState((float) mouseX/kWindowWidth,(float) mouseY/kWindowHeight,lbuttonDown,rbuttonDown);
	for(int i = 0; i < 7; ++i)
	{
		simulation->advance();
	}
	glutPostRedisplay();
}

void drawPoints(std::vector<vec2> points, float r, float g, float b, float a, float size)
{

	glUseProgram (programID);
	GLint loc = glGetUniformLocation(programID, "uColor");
    glUniform4f(loc, r, g, b, a);
	glDisable(GL_POINT_SPRITE);
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

    glPointSize(20.0);

	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    // draw points from the currently bound VAO with current in-use shader
    glDrawArrays (GL_POINTS, 0, points.size());
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//draw to screen
	metaballprocess->renderFrom();
}

void mouse(int button, int state, int x, int y)
{
	mouseX = x;
	mouseY = y;
	if(button == GLUT_RIGHT_BUTTON)
	{
		if(state == GLUT_DOWN)
		{
			std::cout <<"Right button pressed at (" <<x <<","<<y << std::endl;
			rbuttonDown = true;
		}else{
			rbuttonDown = false;
		}

	}
	else if(button == GLUT_LEFT_BUTTON)
	{
		if(state == GLUT_DOWN)
		{
			std::cout <<"Left button pressed at (" <<x <<","<<y << ")"<< std::endl;
			lbuttonDown = true;
		}else{
			lbuttonDown = false;
		}
	}
}

void motion(int x, int y)
{
	mouseX = x;
	mouseY = y;
	if (lbuttonDown){}
		//std::cout << "Mouse dragged with left button at "
		//<< "(" << x << "," << y << ")" << std::endl;
			if(rbuttonDown){}
		//std::cout << "Mouse dragged with right button at "
		//<< "(" << x << "," << y << ")" << std::endl;
}

void keyPressed(unsigned char c, int x, int y)
{
	

	if(c == 'p')
	{
		if(POINTS_MODE)
		{
			POINTS_MODE = false;
		}
		else
		{
			POINTS_MODE = true;
		}
	}

	if(c == 'r')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
	}

	if(c == '1')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '1';
	}

	if(c == '2')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '2';
	}

	if(c == '3')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '3';
	}

	if(c == '4')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '4';
	}

		if(c == '5')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '5';
	}

	if(c == '6')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '6';
	}

	if(c == '7')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '7';
	}

	if(c == '8')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '8';
	}

	if(c == '9')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '9';
	}

	if(c == '0')
	{
		if(POINTS_RELOAD)
		{
			POINTS_RELOAD = false;
		}
		else
		{
			POINTS_RELOAD = true;
		}
		keyPress = '0';
	}
	
}

void valueControllerInit(){
	glClearColor(0, 0, 0, 0);

	glViewport(0, 0, 500, 500);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(0, 500, 0, 500, 1, -1);

	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
}

void valueControllerMouse(int button, int state, int x, int y)
{
	//std::cout << "X: " << x << "      Y: " << (y-300)*-1 << std::endl;

	if(button == GLUT_LEFT_BUTTON && (  x==147||x==148||x==149||x==150||x==151  )){
		//std::cout << "X: " << x << "      Y: " << (y-300)*-1 << std::endl;

		stiffnessValueChanger = (y-300)*-1;
		stiffnessValueChanged = true;
	}

	if(button == GLUT_LEFT_BUTTON && (  x==347||x==348||x==349||x==350||x==351  )){
		//std::cout << "X: " << x << "      Y: " << (y-300)*-1 << std::endl;

		viscosityValueChanger = (y-300)*-1;
		viscosityValueChanged = true;
	}
}

void drawGLString(GLfloat x, GLfloat y, char *textString)
{
    int le;
    int qs;
    
    
    glRasterPos2f(x, y);
    le = (int) strlen(textString);
    for (qs = 0; qs < le; qs++) 
    {
        glutBitmapCharacter(GLUT_BITMAP_8_BY_13, textString[qs]);
        
    }
}