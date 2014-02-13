uniform vec3 uLight;
uniform sampler2D uTexUnit0, uTexUnit1, uTexUnit3;

varying vec3 vNormal;
varying vec3 vPosition;
varying vec2 vTexCoord0;
varying vec2 vTexCoord1;
varying vec2 vTexCoord3;

void main() {
  vec3 tolight = normalize(uLight - vPosition);
  vec3 normal = normalize(vNormal);
  vec4 texColor0 = texture2D(uTexUnit3, vTexCoord3);

  float diffuse = max(0.0, dot(normal, tolight));
  diffuse += max(0.0, dot(normal, tolight));

  gl_FragColor = texColor0 * diffuse;
}
