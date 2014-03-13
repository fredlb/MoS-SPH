#version 400

uniform vec4 uColor;

out vec4 frag_colour;
void main () 
{
	frag_colour = uColor;
}