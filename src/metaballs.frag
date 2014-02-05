
uniform sampler2D fbo_texture;
varying vec2 f_texcoord;

void main(void) {
	vec4 tex = texture2D(fbo_texture, f_texcoord);
	gl_FragColor = vec4(1.0)*step(0.25,tex.r);
}