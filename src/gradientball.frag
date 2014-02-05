#version 400

out vec4 frag_colour;
void main () 
{
	vec2 diff = (gl_PointCoord-0.5);
	float r2 = dot(diff,diff);
	r2 = clamp(r2,0.0,1.0);
	float grad = r2*r2 - r2 + 0.25;
	//float grad = 1/r2;
	frag_colour = vec4(1.0, 0.0, 0.0, grad);
}