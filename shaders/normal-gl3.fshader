#version 130

uniform vec3 uLight;
uniform sampler2D uTexUnit0;
uniform sampler2D uTexUnit1;

in vec3 vNormal;
in vec3 vPosition;
in vec3 vTexCoord0;
in vec3 vTexCoord1;

out vec4 fragColor;

void main() 
{
    vec3 tolight = normalize(uLight - vPosition);
    vec4 texNormal = texture2D(uTexUnit1, vTexCoord1);
    vec3 scaledTexNormal = normalize(2.0 * vec3(texNormal)  - vec3(1., 1., 1.));

    float diffuse = max(0.0, dot(vec3(scaledTexNormal), tolight));
    fragColor = vec4(0.7, 0.7, 0.7, 1) * diffuse;
}