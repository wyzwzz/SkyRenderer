#version 460 core

layout(location = 0) in vec2 iFragTexCoord;
layout(location = 1) in vec2 iScreenCoord;

layout(location = 0) out vec4 oFragColor;

layout(std140,binding = 0) uniform SkyParams{
    vec3 camera_dir;
    float exposure;
    vec3 up;
    int enable_tone_mapping;
    vec3 right;
    float scale;//tan(fov/2)
};

layout(binding = 0) uniform sampler2D SkyView;
uniform float WOverH;

#define PI 3.14159265
#define POSTCOLOR_A 2.51
#define POSTCOLOR_B 0.03
#define POSTCOLOR_C 2.43
#define POSTCOLOR_D 0.59
#define POSTCOLOR_E 0.14

vec3 tonemap(vec3 color)
{
    return (color * (POSTCOLOR_A * color + POSTCOLOR_B))
    / (color * (POSTCOLOR_C * color + POSTCOLOR_D) + POSTCOLOR_E);
}

vec3 toneColorMapping(vec2 seed, in out vec3 color)
{

    color = tonemap(color);
    float rand = fract(sin(dot(seed, vec2(12.9898, 78.233) * 2.0)) * 43758.5453);
    color = 255 * clamp(pow(color, vec3(1 / 2.2)),vec3(0),vec3(1));
    color = all(lessThan(rand.xxx,(color - floor(color)))) ? ceil(color) : floor(color);
    return color / 255.0;

}

vec3 whitePointColorMapping(float exposure,in out vec3 color){
    const vec3 white_point = vec3(1.08241, 0.96756, 0.95003);
    return pow(vec3(1.0) - exp(-color / white_point * exposure), vec3(1.0 / 2.2));
}


void main() {
    //    oFragColor = vec4(iFragTexCoord,1, 1);return;
    float ratio = WOverH;
    vec3 view_dir = normalize(camera_dir + ratio * scale * right * iScreenCoord.x + scale * up * iScreenCoord.y);

    float phi = atan(view_dir.z,view_dir.x);
    float u = (phi < 0 ? phi + 2 * PI : phi) / (2 * PI);


    float theta = asin(view_dir.y);
    float v = 0.5 + 0.5 * sign(theta) * sqrt(abs(theta) / (PI / 2));

    vec3 sky_color = texture(SkyView,vec2(u,v)).rgb;

    sky_color = bool(enable_tone_mapping) ? toneColorMapping(iFragTexCoord,sky_color) : whitePointColorMapping(exposure,sky_color);

    oFragColor = vec4(sky_color,1.0);
}
