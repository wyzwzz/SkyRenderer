#version 460 core

#define PI 3.14159265

layout(location = 0) in vec2 iFragTexCoord;
layout(location = 1) in vec2 iScreenCoord;

layout(location = 0) out vec4 oFragColor;

layout(std140, binding = 0) uniform SkyParams{
    vec3 camera_dir;
    float exposure;
    vec3 up;
    float scale;//tan(fov/2)
    vec3 right;
};

layout(binding = 0) uniform sampler2D SkyView;
layout(binding = 1) uniform sampler2D Cloud;

uniform float WOverH;
uniform int EnableCloud;

vec3 whitePointColorMapping(float exposure, in out vec3 color){
    const vec3 white_point = vec3(1.08241, 0.96756, 0.95003);
    return pow(vec3(1.0) - exp(-color / white_point * exposure), vec3(1.0 / 2.2));
}

void main() {
    float ratio = WOverH;
    vec3 view_dir = normalize(camera_dir + ratio * scale * right * iScreenCoord.x + scale * up * iScreenCoord.y);

    float phi = atan(view_dir.z, view_dir.x);
    float u = (phi < 0 ? phi + 2 * PI : phi) / (2 * PI);

    float theta = asin(view_dir.y);
    float v = 0.5 + 0.5 * sign(theta) * sqrt(abs(theta) / (PI / 2));

    vec3 sky_color = texture(SkyView, vec2(u, v)).rgb;

    sky_color = whitePointColorMapping(exposure, sky_color);

    vec4 cloud_color = vec4(0);

    if (bool(EnableCloud)){
        cloud_color = texture(Cloud, iFragTexCoord);
        cloud_color.rgb *= sky_color;
        cloud_color.rgb = whitePointColorMapping(exposure, cloud_color.rgb);
    }

    oFragColor = vec4(mix(sky_color, cloud_color.rgb, cloud_color.a), 1.0);
}
