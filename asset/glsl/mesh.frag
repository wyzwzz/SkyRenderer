#version 460 core

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

layout(binding = 3) uniform sampler2D AlbedoMap;//optional


vec3 postProcessColor(vec3 color)
{
    vec3 white_point = vec3(1.08241, 0.96756, 0.95003);
    float exposure = 10.0;
    return pow(vec3(1.0) - exp(-color / white_point * exposure), vec3(1.0 / 2.2));
}
vec3 getTransmittance(float h,float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(Transmittance, vec2(u, v)).rgb;
}

void main() {
    //    oFragColor = vec4(iFragNormal,1.f);return;
    vec3 world_pos = iFragPos;
    vec3 normal = normalize(iFragNormal);
    vec2 screen_coord = iScreenCoord.xy / iScreenCoord.w;
    screen_coord = 0.5 + screen_coord * vec2(0.5,-0.5);//ndc differ with image which start at left-up corner

    float z = world_scale * distance(world_pos,view_pos) / max_aeraial_distance;
    vec4 aerial_res = texture(AerialPerspective,vec3(screen_coord,z));
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

    vec3 res = sun_intensity * (
    in_scattering +
    shadow_factor * sun_transmittance * albedo
    * max(0,dot(normal,-sun_dir)) * view_transmittance);

    res = postProcessColor(res);
    oFragColor = vec4(res,1.0);
}
