#include "common.hpp"

struct AtmosphereProperties {
    vec3f rayleigh_scattering = {5.802f, 13.558f, 33.1f};//10^(-6)m^(-1)
    float rayleigh_density_h = 8.f;//km

    float mie_scattering = 3.996f;
    float mie_asymmetry_g = 0.88f;//scalar
    float mie_absorption = 4.4f;
    float mie_density_h = 1.2f;

    vec3f ozone_absorption = {0.65f, 1.881f, 0.085f};//10^(-6)m^(-1)
    // absolute height = ground_radius + ozone_center_h
    float ozone_center_h = 25;//km
    float ozone_width = 30;//km

    float ground_radius = 6360;//km
    float top_atmosphere_radius = 6460;//km
    float padding = 0;

    // to meter
    AtmosphereProperties toStdUnit() const {
        AtmosphereProperties ret = *this;

        ret.rayleigh_scattering *= 1e-6f;
        ret.rayleigh_density_h *= 1e3f;

        ret.mie_scattering *= 1e-6f;
        ret.mie_absorption *= 1e-6f;
        ret.mie_density_h *= 1e3f;

        ret.ozone_absorption *= 1e-6f;
        ret.ozone_center_h *= 1e3f;
        ret.ozone_width *= 1e3f;

        ret.ground_radius *= 1e3f;
        ret.top_atmosphere_radius *= 1e3f;

        return ret;
    }
};

static_assert(sizeof(AtmosphereProperties) == 64, "");


struct alignas(16) CloudProperties {
    vec3 density_to_sigma_s = vec3(1);
    float phase_g = 0.8;
    vec3 density_to_sigma_t = vec3(1);
    int primary_ray_marching_steps = 64;
    vec3 box_min = vec3(-50, 50, -50);
    int secondary_ray_marching_steps = 4;
    vec3 box_max = vec3(50, 60, 50);
    int enable_multi_scattering = 0;//not using...
    float g_c = 0.8;//global coverage control
    float g_d = 1.0;//global density control
    float wc0 = 0.8;
    float wc1 = 0.9;
    float wh = 0.8;
    float shape_tile = 0.03;
    float detail_tile = 0.11;
    float blend_alpha = 0.5;
};

