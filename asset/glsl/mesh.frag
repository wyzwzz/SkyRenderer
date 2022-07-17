#version 460 core
#define PI 3.14159265
layout(location = 0) in vec3 iFragPos;
layout(location = 1) in vec3 iFragNormal;
layout(location = 2) in vec2 iFragTexCoord;
layout(location = 3) in vec4 iScreenCoord;

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
layout(std140,binding = 1) uniform MeshParams{
    vec3 sun_dir;
    float sun_theta;
    vec3 sun_intensity;
    float max_aeraial_distance;
    vec3 view_pos;
    float world_scale;
    mat4 shadow_vp;

};

layout(binding = 0) uniform sampler2D Transmittance;
layout(binding = 1) uniform sampler3D AerialPerspective;
layout(binding = 2) uniform sampler2D ShadowMap;
layout(binding = 3) uniform sampler2D BlueNoise;


layout(binding = 3) uniform sampler2D AlbedoMap;//optional


vec3 postProcessColor(vec3 color)
{
//    return pow(color,vec3(1.0/2.2));
    vec3 white_point = vec3(1.08241, 0.96756, 0.95003);
    float exposure = 10.0;
    return pow(vec3(1.0) - exp(-color / white_point * exposure), vec3(1.0 / 2.2));
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

void main() {
    //    oFragColor = vec4(iFragNormal,1.f);return;
    vec3 world_pos = iFragPos;
    vec3 normal = normalize(iFragNormal);
    vec2 screen_coord = iScreenCoord.xy / iScreenCoord.w;
    screen_coord = 0.5 + screen_coord * vec2(0.5,-0.5);//ndc differ with image which start at left-up corner

    float z = world_scale * distance(world_pos,view_pos) / max_aeraial_distance;
    float noise = texture(BlueNoise,world_pos.xy * 0.01).r;
    vec2 res = vec2(textureSize(AerialPerspective,0).xy);
    vec2 offset = vec2(1) / res * noise * sin(PI * 43786.12f * noise);
    vec4 aerial_res = texture(AerialPerspective,vec3(screen_coord + offset,z));
    //    oFragColor = vec4(screen_coord,z,1.0);return;
    vec3 in_scattering = aerial_res.rgb;
    //    oFragColor = vec4(in_scattering,1.0);return;
    float view_transmittance = aerial_res.w;
    vec3 sun_transmittance = getTransmittance(world_pos.y * world_scale,sun_theta);
    vec3 albedo = vec3(0.03);//texture(AlbedoMap,iFragTexCoord).rgb;

    vec4 shadow_clip = shadow_vp * vec4(world_pos + 0.03 * normal,1);
    vec3 shadow_ndc = shadow_clip.xyz / shadow_clip.w;
    vec2 shadow_uv = 0.5 + 0.5 * shadow_ndc.xy;

    float shadow_factor = 1;
    if(all(equal(clamp(shadow_uv,vec2(0),vec2(1)), shadow_uv))){
        float view_z = shadow_clip.z * 0.5 + 0.5;
        float shadow_z = texture(ShadowMap,shadow_uv).r;
        shadow_factor = float(view_z <= shadow_z);
    }

    vec3 ret = sun_intensity * (
    in_scattering +
    shadow_factor * sun_transmittance * albedo
    * max(0,dot(normal,-sun_dir)) * view_transmittance);

    ret = postProcessColor(ret);
    oFragColor = vec4(ret,1.0);
}
