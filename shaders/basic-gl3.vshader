#version 130

uniform mat4 uProjMatrix;
uniform mat4 uModelViewMatrix;
uniform mat4 uNormalMatrix;

in vec3 aPosition;
in vec3 aNormal;
in vec3 aTangent;

out vec3 vNormal;
out vec3 vPosition;
out vec3 vTangent;

void main() {
  vNormal = vec3(uNormalMatrix * vec4(aNormal, 0.0));

  // send position (eye coordinates) to fragment shader
  vec4 tPosition = uModelViewMatrix * vec4(aPosition, 1.0);
  vPosition = vec3(tPosition);
  vTangent = normalize(vec3(uModelViewMatrix * vec4(vec3(aTangent), 0.0)));
  gl_Position = uProjMatrix * tPosition;
}