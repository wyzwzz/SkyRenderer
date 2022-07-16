#version 460 core
layout(location = 0) in vec3 iVertexPos;

uniform mat4 MVP;
uniform float sun_theta;
uniform float view_height;
layout(std140,binding = 0) uniform AtmosphereProperties{
    vec3 rayleigh_scattering;
    float rayleigh_density_h;

    float mie_scattering;
    float mie_asymmetry_g;
    float mie_absorption;
    float mie_density_h;

    vec3 ozone_absorption;
    float ozone_center_h;
    float ozone_width;

    float ground_radius;
    float top_atmosphere_radius;
};
layout(binding = 0) uniform sampler2D Transmittance;
layout(binding = 1) uniform sampler2D MultiScattering;
layout(location = 0) out vec4 oClipPos;
layout(location = 1) out vec3 oTransmittance;

vec3 getTransmittance(float h,float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(Transmittance, vec2(u, v)).rgb;
}
vec3 getMultiScattering(float h,float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(MultiScattering,vec2(u,v)).rgb;
}
void main() {
    gl_Position = MVP * vec4(iVertexPos,1);
    gl_Position.z = gl_Position.w;
    oClipPos = gl_Position;
    oTransmittance = getTransmittance(view_height,sun_theta);
    oTransmittance += getMultiScattering(view_height,sun_theta);
}
