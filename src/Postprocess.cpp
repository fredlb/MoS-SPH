#include <GL/glew.h>
#include <stdio.h>
#include <iostream>

#include "shader.h"

#include "postprocess.h"

Postprocess::Postprocess(GLsizei screen_width, GLsizei screen_height, const char * fragment_file_path)
{
	/* Framebuffer */
	/* Create back-buffer, used for post-processing */

	/* Texture */
	glActiveTexture(GL_TEXTURE0);
	glGenTextures(1, &fbo_texture);
	glBindTexture(GL_TEXTURE_2D, fbo_texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen_width, screen_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindTexture(GL_TEXTURE_2D, 0);

	/* Depth buffer */
	glGenRenderbuffers(1, &rbo_depth);
	glBindRenderbuffer(GL_RENDERBUFFER, rbo_depth);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, screen_width, screen_height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);

	/* Framebuffer to link everything together */
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, fbo_texture, 0);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo_depth);
	GLenum status;
	if ((status = glCheckFramebufferStatus(GL_FRAMEBUFFER)) != GL_FRAMEBUFFER_COMPLETE) {
		fprintf(stderr, "glCheckFramebufferStatus: error %p", status);
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	/* Vertices */
	GLfloat fbo_vertices[] = {
	-1, -1,
	 1, -1,
	-1,  1,
	 1,  1,
	};
	glGenBuffers(1, &vbo_fbo_vertices);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_fbo_vertices);
	glBufferData(GL_ARRAY_BUFFER, sizeof(fbo_vertices), fbo_vertices, GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	/* Shader program */
	/* Post-processing */
	program_postproc = LoadShader("src/postprocess.vert", "src/metaballs.frag");

	GLchar * attribute_name = "v_coord";
	attribute_v_coord_postproc = glGetAttribLocation(program_postproc, attribute_name);
	if (attribute_v_coord_postproc == -1) {
		fprintf(stderr, "Could not bind attribute %s\n", attribute_name);
	}

	GLchar * uniform_name = "fbo_texture";
	uniform_fbo_texture = glGetUniformLocation(program_postproc, uniform_name);
	if (uniform_fbo_texture == -1) {
		fprintf(stderr, "Could not bind uniform %s\n", uniform_name);
	}
}

Postprocess::~Postprocess()
{
	/* free_resources */
	glDeleteRenderbuffers(1, &rbo_depth);
	glDeleteTextures(1, &fbo_texture);
	glDeleteFramebuffers(1, &fbo);
	glDeleteBuffers(1, &vbo_fbo_vertices);
	glDeleteProgram(program_postproc);
}

void Postprocess::OnReshape(GLsizei screen_width, GLsizei screen_height)
{
	/* onReshape */
	// Rescale FBO and RBO as well
	glBindTexture(GL_TEXTURE_2D, fbo_texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen_width, screen_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glBindTexture(GL_TEXTURE_2D, 0);

	glBindRenderbuffer(GL_RENDERBUFFER, rbo_depth);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT16, screen_width, screen_height);
	glBindRenderbuffer(GL_RENDERBUFFER, 0);
}

void Postprocess::drawTo()
{
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
}

void Postprocess::renderFrom()
{
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glUseProgram(program_postproc);
	glBindTexture(GL_TEXTURE_2D, fbo_texture);
	glUniform1i(uniform_fbo_texture, /*GL_TEXTURE*/0);
	glEnableVertexAttribArray(attribute_v_coord_postproc);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_fbo_vertices);
	glVertexAttribPointer(
		attribute_v_coord_postproc,  // attribute
		2,                  // number of elements per vertex, here (x,y)
		GL_FLOAT,           // the type of each element
		GL_FALSE,           // take our values as-is
		0,                  // no extra data between each position
		0                   // offset of first element
	);
	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	glDisableVertexAttribArray(attribute_v_coord_postproc);
}

GLuint Postprocess::getShaderID()
{
	return program_postproc;
}