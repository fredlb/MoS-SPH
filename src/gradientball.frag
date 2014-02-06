#version 400

out vec4 frag_colour;
void main () 
{
	vec2 diff = (gl_PointCoord-0.5);
	float r2 = dot(diff,diff)*4.0;
	r2 = clamp(r2, 0.0, 0.5);
	float grad = r2*r2 - r2 + 0.25;
	//float grad = (1-r2)*(1-r2);
	//float grad = (1.0 - r2);

	//float grad = 0.05/r2;
	//grad = step(0.2,grad)*grad;

	frag_colour = vec4(grad, 0.0, 0.0, 1.0);
}