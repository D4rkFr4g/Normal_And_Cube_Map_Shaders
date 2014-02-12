uniform vec3 uLight;
uniform sampler2D uTexUnit0;
uniform sampler2D uTexUnit1;

varying vec3 vNormal;
varying vec3 vPosition;
varying vec2 vTexCoord0;
varying vec2 vTexCoord1;

void main() 
{
    vec3 tolight = normalize(uLight - vPosition);
    vec4 texNormal = texture2D(uTexUnit1, vTexCoord1);
    vec3 scaledTexNormal = normalize(2.0 * vec3(texNormal)  - vec3(1., 1., 1.));

    float diffuse = max(0.0, dot(vec3(scaledTexNormal), tolight));
    gl_FragColor = vec4(0.7, 0.7, 0.7, 1.0) * diffuse;
}