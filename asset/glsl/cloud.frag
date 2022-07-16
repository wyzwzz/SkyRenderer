#version 460 core

#define PI 3.14159265

layout(location = 0) in vec2 iUV;
layout(location = 1) in vec2 iScreenCoord;

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
layout(std140,binding = 1) uniform CloudProperties{
    vec3 density_to_sigma_s;
    float phase_g;
    vec3 density_to_sigma_t;
    int primary_ray_marching_steps;
    vec3 box_min;
    int secondary_ray_marching_steps;
    vec3 box_max;
    int enable_multi_scattering;
    float g_c;
    float g_d;
    float wc0;
    float wc1;
    float wh;
    float shape_tile;
    float detail_tile;
    float blend_alpha;
};
layout(std140,binding = 2) uniform CloudParams{
    vec3 camera_dir;
    float scale;//tan(fov/2)
    vec3 up;
    float w_over_h;
    vec3 right;
    vec3 view_pos;
    vec3 sun_dir;
};

layout(binding = 0) uniform sampler2D WeatherMap;
layout(binding = 1) uniform sampler3D ShapeNoise;
layout(binding = 2) uniform sampler3D DetailNoise;
layout(binding = 3) uniform sampler2D BlueNoise;
layout(binding = 4) uniform sampler2D Transmittance;
layout(binding = 5) uniform sampler3D AerialPerspective;


vec3 getTransmittance(float h,float theta){
    float u = h / (top_atmosphere_radius - ground_radius);
    float v = 0.5 + 0.5 * sin(theta);
    return texture(Transmittance, vec2(u, v)).rgb;
}

float remap(float original_value, float original_min, float original_max, float new_min, float new_max)
{
    return new_min + (((original_value - original_min) / (original_max - original_min)) * (new_max - new_min));
}

float saturate(float x){
    return clamp(x,0.0,1.0);
}

vec2 saturate(vec2 x){
    return clamp(x,vec2(0.0),vec2(1.0));
}

vec3 saturate(vec3 x){
    return clamp(x,vec3(0.0),vec3(1.0));
}

bool insideBox(vec3 pos,vec3 box_min,vec3 box_max){
    return pos.x >= box_min.x && pos.x <= box_max.x
    && pos.y >= box_min.y && pos.y <= box_max.y
    && pos.z >= box_min.z && pos.z <= box_max.z;
}

vec2 rayIntersectBox(vec3 box_min,vec3 box_max,vec3 ray_origin,vec3 inv_ray_dir){
    vec3 t0 = (box_min - ray_origin) * inv_ray_dir;
    vec3 t1 = (box_max - ray_origin) * inv_ray_dir;
    vec3 t_min = min(t0,t1);
    vec3 t_max = max(t0,t1);
    float enter_t = max(max(t_min.x,t_min.y),t_min.z);
    float exit_t  = min(min(t_max.x,t_max.y),t_max.z);
    return vec2(enter_t,exit_t);//exit_t > enter_t && exit_t > 0
}

vec3 evalPhaseFunction(float u){

    float g = phase_g;
    float g2 = g * g;
    float u2 = u * u;
    float m = 1 + g2 - 2 * g * u;

    float ph = 3 / (8 * PI) * (1 - g2) * (1 + u2) / ((2 + g2) * m * sqrt(m));

    vec3 pu = vec3(abs(3 / (16 * PI) * (1 + u2)));

    return mix(pu,vec3(ph),blend_alpha) * 0.6 + vec3(0.4);
}


float sampleDensity(vec3 pos){
    vec3 box_size = box_max - box_min;
    vec3 box_center = 0.5 * (box_min + box_max);
    vec2 uv = saturate((pos.xz - box_min.xz) / box_size.xz);

    const float containerEdgeFadeDst = 10;
    float dstFromEdgeX = min(containerEdgeFadeDst, min(pos.x - box_min.x, box_max.x - pos.x));
    float dstFromEdgeZ = min(containerEdgeFadeDst, min(pos.z - box_min.z, box_max.z - pos.z));
    float edgeWeight = min(dstFromEdgeZ, dstFromEdgeX) / containerEdgeFadeDst;

    float wd = texture(WeatherMap,uv).r;

    const float WMc = max(wc0,saturate(g_c - 0.5) * wc1 * 2);
    float ph = saturate((pos.y - box_min.y) / box_size.y);

    float SRb = saturate(remap(ph,0,0.07,0,1));
    float SRt = saturate(remap(ph,wh*0.2,wh,1,0));
    float SA = SRb * SRt;

    float DRb = ph * saturate(remap(ph,0,0.15,0,1));
    float DRt = saturate(remap(ph,0.9,1.0,1,0));
    float DA = g_d * DRb * DRt * wd * 2;

    float SN = texture(ShapeNoise,pos * shape_tile).r;

    SN = remap(SN * SA,1 - g_c * WMc,1, 0, 1) * DA * edgeWeight;

    float DN_fbm = texture(DetailNoise,pos * detail_tile).r;
    float DN_mod = 0.35 * exp(-g_c * 0.75) * mix(DN_fbm,1 - DN_fbm,saturate(ph * 5));

    float d = saturate(remap(SN,DN_mod,1,0,1));
    return d;
}

vec3 singleScattering(vec3 pos){
    //ray marching toward sun
    vec3 ray_dir = -sun_dir;
    vec2 intersect_t = rayIntersectBox(box_min,box_max,pos,1.0 / (ray_dir ));
    float ray_march_dist = max(0,intersect_t.y);
    float dt = ray_march_dist / secondary_ray_marching_steps;
    vec3 sum_sigma_t = vec3(0);
    for(int i = 0; i < secondary_ray_marching_steps; ++i){
        vec3 ith_pos = pos + (i + 0.5) * ray_dir * dt;
        sum_sigma_t += sampleDensity(ith_pos) * density_to_sigma_t * dt;
    }
    //todo pos correct with world scale
    return getTransmittance(pos.y,asin(ray_dir.y)) * exp(-sum_sigma_t);
}

vec4 cloudRayMarching(vec3 start_pos,vec3 ray_dir,float max_ray_advance_dist){

    vec3 sum_sigma_t = vec3(0);
    float dt = max_ray_advance_dist / primary_ray_marching_steps;
    vec3 accu_radiance = vec3(0);

    float blue_noise = texture(BlueNoise,iUV).r;
    float advanced_dist = dt * blue_noise ;
    for(int i = 0; i < primary_ray_marching_steps; ++i){
        vec3 ith_pos = start_pos + advanced_dist * ray_dir;

        vec3 ith_single_scattering = singleScattering(ith_pos);

        float density = sampleDensity(ith_pos);

        sum_sigma_t +=  density * density_to_sigma_t * dt;

        vec3 sigma_s = density *  density_to_sigma_s;

        vec3 rho = evalPhaseFunction(dot(sun_dir,-ray_dir));

        vec3 ith_transmittance = vec3(exp(-sum_sigma_t));

        accu_radiance += ith_single_scattering * sigma_s * rho *  ith_transmittance * dt;

        if( all( lessThan( ith_transmittance,vec3(0.01) ) ) )
            break;

        advanced_dist += dt;
    }

    return vec4(accu_radiance, 1 - exp(-sum_sigma_t));
}

void main() {

    vec3 view_dir = normalize(camera_dir + w_over_h * scale * right * iScreenCoord.x + scale * up * iScreenCoord.y);

    vec2 intersect_t = rayIntersectBox(box_min,box_max,view_pos,1.0 / view_dir);


    //todo use aerial lut

    vec4 cloud_color = vec4(0);
    if(intersect_t.y > 0 && intersect_t.y > intersect_t.x){
        cloud_color = cloudRayMarching(view_pos + view_dir * max(0,intersect_t.x),view_dir,intersect_t.y - max(0,intersect_t.x));
    }

    oFragColor = cloud_color;
}
