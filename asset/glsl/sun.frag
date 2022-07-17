#version 460 core

layout(location = 0) in vec4 iClipPos;
layout(location = 1) in vec3 iTransmittance;

layout(location = 0) out vec4 oFragColor;


uniform vec3 SunRadiance;

vec3 postProcessColor(vec2 seed, vec3 color)
{
    vec3 white_point = vec3(1.08241, 0.96756, 0.95003);
    float exposure = 12.0;
    return pow(vec3(1.0) - exp(-color / white_point * exposure), vec3(1.0 / 2.2));

}

void main() {
    vec3 color;
    color = iTransmittance * SunRadiance;
    color = postProcessColor(iClipPos.xy/iClipPos.w, color);
    oFragColor = vec4(color, 1.0);
}
