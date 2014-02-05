
uniform sampler2D fbo_texture;
varying vec2 f_texcoord;

void main(void) {
	vec4 tex = texture2D(fbo_texture, f_texcoord);
	gl_FragColor = vec4(step(0.2,tex.r),step(0.4,tex.r),step(0.6,tex.r),1.0);
}