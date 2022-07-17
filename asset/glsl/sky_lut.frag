#version 460 core

layout(location = 0) in vec2 iFragTexCoord;

layout(location = 0) out vec4 oFragColor;

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

layout(std140,binding = 1) uniform SkyParams{
    vec3 view_pos;
    int ray_march_step_count;

    vec3 sun_dir;//from sun to ...
    int enable_multiscattering;

    vec3 sun_intensity;
};

layout(binding = 0) uniform sampler2D Transmittance;
layout(binding = 1) uniform sampler2D MultiScattering;

#define PI 3.14159265

bool findClosestIntersectionWithCircle(vec2 o, vec2 d, float R, out float t){
    float A = dot(d, d);
    float B = 2 * dot(o, d);
    float C = dot(o, o) - R * R;
    float delta = B * B - 4 * A * C;
    if(delta < 0)
    return false;
    t = (-B + (C <= 0 ? sqrt(delta) : -sqrt(delta))) / (2 * A);
    return (C <= 0) || (B <= 0);
}
void getSigmaST(float h,out vec3 sigma_s,out vec3 sigma_t){
    vec3 rayleigh = rayleigh_scattering * exp(-h / rayleigh_density_h);
    vec3 mie_s = vec3(mie_scattering) * exp(-h / mie_density_h);
    vec3 mie_t = vec3(mie_scattering + mie_absorption) * exp(-h / mie_density_h);
    //d^{ozone}(h) = max(0,1 - frac{|h - 25|}{15})
    vec3 ozone = ozone_absorption * max(0.0,1 - abs(h - ozone_center_h) / (ozone_width * 0.5));
    sigma_s = rayleigh + mie_s;
    sigma_t = rayleigh + mie_t + ozone;
}
bool hasIntersectionWithSphere(vec3 o, vec3 d, float R){
    float A = dot(d, d);
    float B = 2 * dot(o, d);
    float C = dot(o, o) - R * R;
    float delta = B * B - 4 * A * C;
    return (delta >= 0) && ((C <= 0) || (B <= 0));
}
vec3 evalPhaseFunction(float h,float u){
    vec3 sRayleigh = rayleigh_scattering * exp(-h / rayleigh_density_h);
    float sMie = mie_scattering * exp(-h / mie_density_h);
    vec3 s = sRayleigh + sMie;

    float g = mie_asymmetry_g, g2 = g * g, u2 = u * u;
    float pRayleigh = 3 / (16 * PI) * (1 + u2);

    float m = 1 + g2 - 2 * g * u;
    float pMie = 3 / (8 * PI) * (1 - g2) * (1 + u2) / ((2 + g2) * m * sqrt(m));

    vec3 result;
    result.x = s.x > 0 ? (pRayleigh * sRayleigh.x + pMie * sMie) / s.x : 0;
    result.y = s.y > 0 ? (pRayleigh * sRayleigh.y + pMie * sMie) / s.y : 0;
    result.z = s.z > 0 ? (pRayleigh * sRayleigh.z + pMie * sMie) / s.z : 0;
    return result;
}
//----------------------------------------
float safeSqrt(float a) {
    return sqrt(max(a, 0.0));
}

float clampCosine(float mu) {
    return clamp(mu, -1.0, 1.0);
}

float clampDistance(float d) {
    return max(d, 0.0);
}

float clampRadius(float r) {
    return clamp(r, ground_radius, top_atmosphere_radius);
}
//余弦定理
float distanceToTopAtmosphereBoundary(float r,float mu){
    float discriminant = r * r * (mu * mu - 1.0) + top_atmosphere_radius * top_atmosphere_radius;
    return clampDistance(-r * mu + safeSqrt(discriminant));
}
float getTextureCoordFromUnitRange(float x, int texture_size) {
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}
vec2 getTransmittanceTextureUvFromRMu(float r, float mu,in ivec2 res){
    float H =  sqrt(top_atmosphere_radius * top_atmosphere_radius - ground_radius * ground_radius);
    float rho = safeSqrt(r * r - ground_radius * ground_radius);
    float d = distanceToTopAtmosphereBoundary(r,mu);
    float d_min = (top_atmosphere_radius - r);
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return vec2(getTextureCoordFromUnitRange(x_mu,res.x),getTextureCoordFromUnitRange(x_r,res.y));
}

// theta is view with horizon
vec3 getTransmittance(float h,float theta){
    //    float u = h / (top_atmosphere_radius - ground_radius);
    //    float v = 0.5 + 0.5 * sin(theta);
    float r = h * 0.99  + ground_radius;
    float mu = cos(PI / 2 - theta);
    vec2 uv = getTransmittanceTextureUvFromRMu(r,mu,textureSize(Transmittance,0));
    return texture(Transmittance, uv).rgb;
}
vec3 getMultiScattering(float h,float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(MultiScattering,vec2(u,v)).rgb;
}
void main() {
    //non-linear mapping
    //    oFragColor = vec4(iFragTexCoord,0.5, 1);return;
    float phi = 2 * PI * iFragTexCoord.x;
    float vm = 2 * iFragTexCoord.y - 1;
    float theta = sign(vm) * (PI / 2) * vm * vm;
    float sin_theta = sin(theta);
    float cos_theta = cos(theta);

    vec3 view_dir = vec3(cos(phi) * cos_theta,sin_theta,sin(phi) * cos_theta);

    vec2 planet_pos = vec2(0,view_pos.y + ground_radius);
    vec2 planet_dir = vec2(cos_theta,sin_theta);

    float end_t = 0;
    if(!findClosestIntersectionWithCircle(planet_pos,planet_dir,ground_radius,end_t)){
        findClosestIntersectionWithCircle(planet_pos,planet_dir,top_atmosphere_radius,end_t);
    }

    float phase_u = dot(-sun_dir,view_dir);

    vec3 in_scattering_sum = vec3(0);
    vec3 sigam_t_sum = vec3(0);
    float t = 0;
    float dt = end_t / ray_march_step_count;
    float half_dt = 0.5 * dt;
    for(int i = 0; i < ray_march_step_count; ++i){
        float mid_t = t + half_dt;
        t += dt;
        vec3 ith_pos = vec3(0,view_pos.y + ground_radius,0) + view_dir * mid_t;
        float ith_h = length(ith_pos) - ground_radius;

        vec3 ith_sigma_s,ith_sigma_t;
        getSigmaST(ith_h,ith_sigma_s,ith_sigma_t);

        vec3 ith_transmittance = exp(-(sigam_t_sum + ith_sigma_t * half_dt));

        float sun_theta = PI / 2 - acos(dot(-sun_dir,normalize(ith_pos)));
        if(!hasIntersectionWithSphere(ith_pos,-sun_dir,ground_radius)){
            vec3 rho = evalPhaseFunction(ith_h,phase_u);
            vec3 ith_sun_transmittance = getTransmittance(ith_h,sun_theta);
            in_scattering_sum += ith_sun_transmittance * rho * ith_sigma_s * ith_transmittance * dt;
        }

        if(bool(enable_multiscattering)){
            in_scattering_sum += getMultiScattering(ith_h,sun_theta) * ith_sigma_s * ith_transmittance * dt;
        }
        sigam_t_sum += ith_sigma_t * dt;
    }
    oFragColor = vec4(in_scattering_sum * sun_intensity, 1);
}
