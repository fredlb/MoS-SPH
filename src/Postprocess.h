//#pragma once
#include <GL/glew.h>

class Postprocess
{
public:
	Postprocess(GLsizei screen_width, GLsizei screen_height, const char * fragment_file_path);
	~Postprocess();
	void OnReshape(GLsizei screen_width, GLsizei screen_height);
	void drawTo();
	void renderFrom();
	GLuint getShaderID();
private:
	GLuint fbo, fbo_texture, rbo_depth; //framebuffer
	GLuint vbo_fbo_vertices; //vertices
	GLuint program_postproc, attribute_v_coord_postproc, uniform_fbo_texture; //program
};