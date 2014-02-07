uniform sampler2D fbo_texture;
varying vec2 f_texcoord;

const vec2 resolution = vec2(800.0, 600.0);

void main(void) {
	float height = texture2D(fbo_texture, f_texcoord).r; //height at fragment
	height = step(0.2,height)*height;
	

	vec2 uv = f_texcoord/resolution;

	vec3 h = vec3(1.0 / resolution, 0.0);
	float hx1 = texture2D(fbo_texture, f_texcoord+h.xz).r;
	float hy1 = texture2D(fbo_texture, f_texcoord+h.zy).r;

	vec3 delta;
	delta.x = (hx1 - height)*resolution.x;
	delta.y = (hy1 - height)*resolution.y;
	delta.z = 40.0;

	vec3 normal = normalize(delta);
	
	// light direction
	vec3 light = normalize(vec3(-1.0, -1.0, 0.5));
	
	// smooth border
	float alpha = smoothstep(0.0, 0.3, height)-0.4;
	
	// reflected light vector
	vec3 ref = reflect(normalize(vec3(uv, -1.0)), normal);
	
	// diffuse light term
	float diff = dot(normal, light) * 0.6 + 0.4;
	
	// combine colors
	vec3 col = vec3(0.0, 0.75, 0.95);
	col = col * diff;
	
	//gl_FragColor = vec4(normal,step(0.01, height));
	gl_FragColor = vec4(col,alpha);
	//gl_FragColor = vec4(1.0)*height;
}