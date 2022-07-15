#version 460 core

layout(location = 0) in vec3 iFragPos;
layout(location = 1) in vec3 iFragNormal;
layout(location = 2) in vec2 iFragTexCoord;
layout(location = 3) in vec4 iScreenCoord;

layout(location = 0) out vec4 oFragColor;

layout(binding = 0) uniform sampler2D Transmittance;
layout(binding = 1) uniform sampler3D AerialPerspective;
layout(binding = 2) uniform sampler2D ShadowMap;

void main(){
    oFragColor = vec4(iFragNormal * 0.5 + 0.5,1.0);
}