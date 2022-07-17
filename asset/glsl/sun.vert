#version 460 core

#define PI 3.14159265

layout(location = 0) in vec3 iVertexPos;

layout(location = 0) out vec4 oClipPos;
layout(location = 1) out vec3 oTransmittance;

layout(std140, binding = 0) uniform AtmosphereProperties{
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

uniform mat4 MVP;
uniform float sun_theta;
uniform float view_height;

float safeSqrt(float a) {
    return sqrt(max(a, 0.0));
}

float clampDistance(float d) {
    return max(d, 0.0);
}

float distanceToTopAtmosphereBoundary(float r, float mu){
    float discriminant = r * r * (mu * mu - 1.0) + top_atmosphere_radius * top_atmosphere_radius;
    return clampDistance(-r * mu + safeSqrt(discriminant));
}

float getTextureCoordFromUnitRange(float x, int texture_size) {
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

vec2 getTransmittanceTextureUvFromRMu(float r, float mu, in ivec2 res){
    float H =  sqrt(top_atmosphere_radius * top_atmosphere_radius - ground_radius * ground_radius);
    float rho = safeSqrt(r * r - ground_radius * ground_radius);
    float d = distanceToTopAtmosphereBoundary(r, mu);
    float d_min = (top_atmosphere_radius - r);
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return vec2(getTextureCoordFromUnitRange(x_mu, res.x), getTextureCoordFromUnitRange(x_r, res.y));
}

// theta is view with horizon
vec3 getTransmittance(float h, float theta){
    float r = h * 0.99  + ground_radius;
    float mu = cos(PI / 2 - theta);
    vec2 uv = getTransmittanceTextureUvFromRMu(r, mu, textureSize(Transmittance, 0));
    return texture(Transmittance, uv).rgb;
}

vec3 getMultiScattering(float h, float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(MultiScattering, vec2(u, v)).rgb;
}

bool hasIntersectionWithSphere(vec3 o, vec3 d, float R){
    float A = dot(d, d);
    float B = 2 * dot(o, d);
    float C = dot(o, o) - R * R;
    float delta = B * B - 4 * A * C;
    return (delta >= 0) && ((C <= 0) || (B <= 0));
}

void main() {
    gl_Position = MVP * vec4(iVertexPos, 1);
    gl_Position.z = gl_Position.w;
    oClipPos = gl_Position;
    vec3 p = vec3(0, view_height + ground_radius, 0);
    vec3 d = vec3(cos(sun_theta), sin(sun_theta), 0);
    if (hasIntersectionWithSphere(p, d, ground_radius))
        oTransmittance = vec3(0);
    else
        oTransmittance = getTransmittance(view_height, sun_theta);

    oTransmittance += getMultiScattering(view_height, sun_theta);
}
