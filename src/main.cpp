#include "sky.hpp"

#ifdef NDEBUG
#define set_uniform_var set_uniform_var_unchecked
#endif

inline vec3i getGroupSize(int x, int y = 1, int z = 1) {
    constexpr int group_thread_size_x = 16;
    constexpr int group_thread_size_y = 16;
    constexpr int group_thread_size_z = 16;
    const int group_size_x = (x + group_thread_size_x - 1) / group_thread_size_x;
    const int group_size_y = (y + group_thread_size_y - 1) / group_thread_size_y;
    const int group_size_z = (z + group_thread_size_z - 1) / group_thread_size_z;
    return {group_size_x, group_size_y, group_size_z};
}

class TransmittanceGenerator {
public:
    void initialize() {
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/transmittance.comp"));
        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);
    }

    void generate(const AtmosphereProperties &ap, const vec2i &lut_size) {
        atmosphere_buffer.set_buffer_data(&ap);
        atmosphere_buffer.bind(0);

        lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1, GL_RGBA32F, lut_size.x, lut_size.y);

        lut.bind_image(0, 0, GL_READ_WRITE, GL_RGBA32F);
        shader.bind();

        auto group_size = getGroupSize(lut_size.x, lut_size.y);
        GL_EXPR(glDispatchCompute(group_size.x, group_size.y, 1));

        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));

        shader.unbind();
    }

    const texture2d_t &getLUT() const {
        return lut;
    }

private:
    program_t shader;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    texture2d_t lut;
};

class MultiScatteringGenerator {
public:
    void initialize() {
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/multiscattering.comp"));
        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(&params, GL_STATIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    }

    void setSample(int ray_march_steps, int dir_samples) {
        params.ray_march_step_count = ray_march_steps;
        if (params.dir_sample_count != dir_samples) {
            params.dir_sample_count = dir_samples;
            dir_samples_buffer.destroy();
            dir_samples_buffer.initialize_handle();
            auto dir_samples_arr = getPoissonDiskSamples(dir_samples);
            dir_samples_buffer.initialize_buffer_data(dir_samples_arr.data(), dir_samples_arr.size() * sizeof(vec2f),
                                                      GL_DYNAMIC_STORAGE_BIT);
        }
    }

    void generate(const AtmosphereProperties &ap, const vec2i &lut_size,
                  const vec3f &albedo,
                  const texture2d_t &transmittance) {
        atmosphere_buffer.set_buffer_data(&ap);
        atmosphere_buffer.bind(0);

        lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1, GL_RGBA32F, lut_size.x, lut_size.y);
        lut.bind_image(0, 0, GL_READ_WRITE, GL_RGBA32F);

        params.ground_albedo = albedo;
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(1);

        dir_samples_buffer.bind(0);

        transmittance.bind(0);
        linear_sampler.bind(0);

        shader.bind();

        auto group_size = getGroupSize(lut_size.x, lut_size.y);
        GL_EXPR(glDispatchCompute(group_size.x, group_size.y, 1));

        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));

        shader.unbind();
    }

    const texture2d_t &getLUT() const {
        return lut;
    }

private:
    struct Params {
        vec3f ground_albedo = vec3(0.3f);
        int dir_sample_count = 0;
        vec3f sun_intensity = vec3f(1);//todo
        int ray_march_step_count = 0;
    } params;
    std140_uniform_block_buffer_t<Params> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    storage_buffer_t<vec2f> dir_samples_buffer;

    program_t shader;
    texture2d_t lut;

    sampler_t linear_sampler;
};

class ShadowGenerator {
public:
    void initialize() {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/mesh.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/depth.frag"));
        shadow.fbo.initialize_handle();
        shadow.rbo.initialize_handle();
        shadow.rbo.set_format(GL_DEPTH32F_STENCIL8, 4096, 4096);
        shadow.fbo.attach(GL_DEPTH_STENCIL_ATTACHMENT, shadow.rbo);
        shadow.tex.initialize_handle();
        shadow.tex.initialize_texture(1, GL_R32F, 4096, 4096);
        shadow.fbo.attach(GL_COLOR_ATTACHMENT0, shadow.tex);
        assert(shadow.fbo.is_complete());
    }

    void begin() {
        shadow.fbo.bind();
        GL_EXPR(glViewport(0, 0, 4096, 4096));
        GL_EXPR(glDrawBuffer(GL_COLOR_ATTACHMENT0));
        shadow.fbo.clear_color_depth_buffer();
        shader.bind();
    }

    void generate(const vertex_array_t &vao, GLsizei count, const mat4 &model, const mat4 &proj_view) {
        shader.set_uniform_var("Model", model);
        shader.set_uniform_var("ProjView", proj_view);
        vao.bind();
        GL_EXPR(glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, nullptr));
        vao.unbind();
    }

    void end() {
        shader.unbind();
        shadow.fbo.unbind();
    }

    const texture2d_t &getShadowMap() const {
        return shadow.tex;
    }

private:
    program_t shader;
    struct {
        framebuffer_t fbo;
        renderbuffer_t rbo;
        texture2d_t tex;
    } shadow;
};

class SkyLUTGenerator {
public:
    void initialize(const vec2i &res) {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/quad.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/sky_lut.frag"));

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        sky.fbo.initialize_handle();
        resize(res);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr, GL_DYNAMIC_DRAW);

        vao.initialize_handle();
    }

    void resize(const vec2i &res) {
        sky.rbo.destroy();
        sky.rbo.initialize_handle();
        sky.rbo.set_format(GL_DEPTH24_STENCIL8, res.x, res.y);
        sky.fbo.attach(GL_DEPTH_STENCIL_ATTACHMENT, sky.rbo);
        sky.lut.destroy();
        sky.lut.initialize_handle();
        sky.lut.initialize_texture(1, GL_RGBA32F, res.x, res.y);
        sky.fbo.attach(GL_COLOR_ATTACHMENT0, sky.lut);
        assert(sky.fbo.is_complete());
        sky_res = res;
    }

    void setAtmosphere(const AtmosphereProperties &atmos) {
        atmosphere_buffer.set_buffer_data(&atmos);
    }

    void setSun(const vec3f &sun_dir, const vec3f &sun_radiance) {
        params.sun_dir = sun_dir;
        params.sun_radiance = sun_radiance;
    }

    void generate(const vec3f &view_pos, int ray_march_steps, bool enable_multi_scattering,
                  const texture2d_t &transmittance,
                  const texture2d_t &multi_scattering) {
        params.view_pos = view_pos;
        params.ray_march_step_count = ray_march_steps;
        params.enable_multi_scattering = static_cast<int>(enable_multi_scattering);
        params_buffer.set_buffer_data(&params);
        atmosphere_buffer.bind(0);
        params_buffer.bind(1);

        transmittance.bind(0);
        multi_scattering.bind(1);
        linear_sampler.bind(0);
        linear_sampler.bind(1);

        sky.fbo.bind();
        GL_EXPR(glViewport(0, 0, sky_res.x, sky_res.y));
        GL_EXPR(glDrawBuffer(GL_COLOR_ATTACHMENT0));
        sky.fbo.clear_color_depth_buffer();
        shader.bind();
        vao.bind();

        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawArrays(GL_TRIANGLE_STRIP, 0, 4));
        GL_EXPR(glDepthFunc(GL_LESS));

        vao.unbind();
        shader.unbind();
        sky.fbo.unbind();
    }

    const texture2d_t &getLUT() const {
        return sky.lut;
    }

private:
    struct Params {
        vec3f view_pos;
        int ray_march_step_count;
        vec3f sun_dir;
        int enable_multi_scattering;
        vec3 sun_radiance;
        float pad = 0.f;
    } params;
    std140_uniform_block_buffer_t<Params> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;

    program_t shader;
    sampler_t linear_sampler;

    struct {
        framebuffer_t fbo;
        renderbuffer_t rbo;
        texture2d_t lut;
    } sky;
    vertex_array_t vao;
    vec2i sky_res;
};

class AerialLUTGenerator {
public:
    void initialize(const vec3i &res) {
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/aerial_lut.comp"));

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr, GL_DYNAMIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        nearest_sampler.initialize_handle();
        nearest_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        resize(res);
    }

    void resize(const vec3i &res) {
        lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1, GL_RGBA32F, res.x, res.y, res.z);
        lut_size = res;
        params.slice_z_count = res.z;
    }

    void setAtmosphere(const AtmosphereProperties &atmos) {
        atmosphere_buffer.set_buffer_data(&atmos);
    }

    void setCamera(const mat4 &inv_proj_view) {
        const vec4f A0 = inv_proj_view * vec4f(-1, 1, 0.2f, 1);
        const vec4f A1 = inv_proj_view * vec4f(-1, 1, 0.5f, 1);

        const vec4f B0 = inv_proj_view * vec4f(1, 1, 0.2f, 1);
        const vec4f B1 = inv_proj_view * vec4f(1, 1, 0.5f, 1);

        const vec4f C0 = inv_proj_view * vec4f(-1, -1, 0.2f, 1);
        const vec4f C1 = inv_proj_view * vec4f(-1, -1, 0.5f, 1);

        const vec4f D0 = inv_proj_view * vec4f(1, -1, 0.2f, 1);
        const vec4f D1 = inv_proj_view * vec4f(1, -1, 0.5f, 1);

        params.frustum_a = (A1 / A1.w - A0 / A0.w).xyz();
        params.frustum_b = (B1 / B1.w - B0 / B0.w).xyz();
        params.frustum_c = (C1 / C1.w - C0 / C0.w).xyz();
        params.frustum_d = (D1 / D1.w - D0 / D0.w).xyz();
    }

    void setSun(const vec3 &sun_dir, const mat4 &sun_proj_view) {
        params.sun_dir = sun_dir;
        params.sun_theta = std::asin(-sun_dir.y);
        params.shadow_vp = sun_proj_view;
    }

    void setAerialPerspective(float world_scale,
                              float max_aerial_distance,
                              int ray_march_steps_per_slice,
                              bool enable_shadow,
                              bool enable_multi_scattering) {
        params.world_scale = world_scale;
        params.max_view_distance = max_aerial_distance;
        params.ray_march_step_count_per_slice = ray_march_steps_per_slice;
        params.enable_shadow = static_cast<int>(enable_shadow);
        params.enable_multiscattering = static_cast<int>(enable_multi_scattering);
    }

    void generate(
            const texture2d_t &transmittance,
            const texture2d_t &multi_scattering,
            const texture2d_t &shadow_map,
            const vec3f &view_pos) {
        atmosphere_buffer.bind(0);

        params.view_pos = view_pos;
        params.view_height = view_pos.y * params.world_scale;
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(1);

        transmittance.bind(0);
        multi_scattering.bind(1);
        shadow_map.bind(2);
        linear_sampler.bind(0);
        linear_sampler.bind(1);
        nearest_sampler.bind(2);

        lut.bind_image(0, 0, GL_READ_WRITE, GL_RGBA32F);

        shader.bind();
        auto group_size = getGroupSize(lut_size.x, lut_size.y, 1);
        GL_EXPR(glDispatchCompute(group_size.x, group_size.y, 1));
        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));
        shader.unbind();

    }

    const texture3d_t &getLUT() const {
        return lut;
    }

private:
    struct AerialParams {
        vec3f sun_dir;
        float sun_theta;
        vec3f frustum_a;
        float max_view_distance;
        vec3f frustum_b;
        int ray_march_step_count_per_slice;
        vec3f frustum_c;
        int enable_multiscattering;
        vec3f frustum_d;
        float view_height;
        vec3f view_pos;
        int enable_shadow;
        mat4 shadow_vp;
        float world_scale;
        int slice_z_count = 32;
        float padding[2];
    } params;
    std140_uniform_block_buffer_t<AerialParams> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;

    texture3d_t lut;
    vec3i lut_size;

    sampler_t linear_sampler;
    sampler_t nearest_sampler;

    program_t shader;
};

class MeshRenderer {
public:
    void initialize(const vec2i &res) {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/mesh.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/mesh.frag"));

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr, GL_DYNAMIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
        nearest_sampler.initialize_handle();
        nearest_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_S, GL_REPEAT);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_T, GL_REPEAT);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_R, GL_REPEAT);

        blue_noise.initialize_handle();
        blue_noise.initialize_format_and_data(1, GL_RGBA8,
                                              wzz::image::load_rgba_from_file("asset/bluenoise.png"));

        resize(res);
    };

    void resize(const vec2i &res) {
        this->res = res;
    }

    void setAtmosphere(const AtmosphereProperties &atmos) {
        atmosphere_buffer.set_buffer_data(&atmos);
    }

    void setSun(const vec3f &sun_dir, const vec3f &sun_radiance, const mat4 &sun_proj_view) {
        params.sun_dir = sun_dir;
        params.sun_theta = std::asin(-sun_dir.y);
        params.shadow_vp = sun_proj_view;
        params.sun_intensity = sun_radiance;
    }

    void setAerialPerspective(float max_aerial_dist, float world_scale) {
        params.max_aerial_distance = max_aerial_dist;
        params.world_scale = world_scale;
    }

    void begin(const texture2d_t &transmittance,
               const texture3d_t &aerial_lut,
               const texture2d_t &shadow_map,
               const mat4 &camera_proj_view) {
        transmittance.bind(0);
        aerial_lut.bind(1);
        shadow_map.bind(2);
        blue_noise.bind(3);
        linear_sampler.bind(0);
        linear_sampler.bind(1);
        nearest_sampler.bind(2);
        nearest_sampler.bind(3);

        GL_EXPR(glViewport(0, 0, res.x, res.y));
        shader.bind();
        shader.set_uniform_var("ProjView", camera_proj_view);
    }

    void render(const mat4 &model,
                const vertex_array_t &vao,
                int draw_count,
                const vec3f &view_pos) {
        shader.set_uniform_var("Model", model);

        atmosphere_buffer.bind(0);

        params.view_pos = view_pos;
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(1);

        vao.bind();
        GL_EXPR(glDrawElements(GL_TRIANGLES, draw_count, GL_UNSIGNED_INT, nullptr));
        vao.unbind();
    }

    void end() {
        shader.unbind();
    }

private:
    program_t shader;
    struct MeshParams {
        vec3f sun_dir;
        float sun_theta;
        vec3f sun_intensity;
        float max_aerial_distance;
        vec3f view_pos;
        float world_scale;
        mat4 shadow_vp;
    } params;
    std140_uniform_block_buffer_t<MeshParams> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;

    vec2i res;

    sampler_t linear_sampler;
    sampler_t nearest_sampler;

    texture2d_t blue_noise;
};

class SkyViewRenderer {
public:
    void initialize(const vec2i &res) {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/quad.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/sky.frag"));

        vao.initialize_handle();

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr, GL_DYNAMIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        nearest_sampler.initialize_handle();
        nearest_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        resize(res);
    }

    void resize(const vec2i &res) {
        view_res = res;
        shader.bind();
        shader.set_uniform_var("WOverH", res.x * 1.f / res.y);
        shader.unbind();
    }

    void enableCloud(bool enable_cloud) {
        shader.bind();
        shader.set_uniform_var("EnableCloud", static_cast<int>(enable_cloud));
        shader.unbind();
    }

    void render(const texture2d_t &sky_lut,
                const texture2d_t &cloud,
                const vec3 &camera_dir,
                const vec3 &camera_right,
                float camera_fov_rad,
                float exposure) {
        params.camera_dir = camera_dir;
        params.right = camera_right;
        params.up = wzz::math::cross(camera_right, camera_dir).normalized();
        params.exposure = exposure;
        params.scale = std::tan(0.5f * camera_fov_rad);
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(0);

        GL_EXPR(glViewport(0, 0, view_res.x, view_res.y));
        shader.bind();
        sky_lut.bind(0);
        cloud.bind(1);
        linear_sampler.bind(0);
        linear_sampler.bind(1);

        vao.bind();

        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawArrays(GL_TRIANGLE_STRIP, 0, 4));
        GL_EXPR(glDepthFunc(GL_LESS));

        vao.unbind();
        shader.unbind();
    }

private:
    struct SkyParams {
        vec3f camera_dir;
        float exposure;
        vec3 up;
        float scale;//tan(fov/2)
        vec3 right;
        float pad;
    } params;
    std140_uniform_block_buffer_t<SkyParams> params_buffer;
    program_t shader;
    vertex_array_t vao;
    vec2i view_res;
    sampler_t linear_sampler;
    sampler_t nearest_sampler;
};

class SunDiskRenderer {
public:
    void initialize(const vec2i &res) {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/sun.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/sun.frag"));

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        auto sphere = MakeSphere();
        sun.vao.initialize_handle();
        sun.pos.initialize_handle();
        sun.ebo.initialize_handle();
        sun.pos.reinitialize_buffer_data(sphere.positions.data(), sphere.positions.size(), GL_STATIC_DRAW);
        sun.vao.bind_vertex_buffer_to_attrib(attrib_var_t<vec3f>(0), sun.pos, 0);
        sun.vao.enable_attrib(attrib_var_t<vec3f>(0));
        sun.ebo.reinitialize_buffer_data(sphere.indices.data(), sphere.indices.size(), GL_STATIC_DRAW);
        sun.vao.bind_index_buffer(sun.ebo);

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        resize(res);
    }

    void resize(const vec2i &res) {
        this->res = res;
    }

    void setAtmosphere(const AtmosphereProperties &atmos) {
        atmosphere_buffer.set_buffer_data(&atmos);
    }

    void setSun(const vec3f &sun_dir, const vec3f &sun_radiance) {
        shader.bind();
        shader.set_uniform_var("sun_theta", std::asin(-sun_dir.y));
        shader.set_uniform_var("SunRadiance", sun_radiance);
        shader.unbind();
    }

    void render(const texture2d_t &transmittance,
                const texture2d_t &multi_scattering,
                const mat4 &mvp,
                float view_height) {
        shader.bind();
        sun.vao.bind();
        shader.set_uniform_var("MVP", mvp);
        shader.set_uniform_var("view_height", view_height);
        transmittance.bind(0);
        multi_scattering.bind(1);
        linear_sampler.bind(0);
        linear_sampler.bind(1);
        atmosphere_buffer.bind(0);
        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawElements(GL_TRIANGLE_STRIP, sun.ebo.index_count(), GL_UNSIGNED_INT, nullptr));
        GL_EXPR(glDepthFunc(GL_LESS));
        sun.vao.unbind();
        shader.unbind();
    }

private:
    program_t shader;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    sampler_t linear_sampler;
    vec2i res;
    struct {
        vertex_array_t vao;
        vertex_buffer_t<vec3f> pos;
        index_buffer_t<uint32_t> ebo;
    } sun;
};

class CloudRenderer {
public:
    void initialize(const vec2i &res) {
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/quad.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/cloud.frag"));

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S, GL_REPEAT);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T, GL_REPEAT);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R, GL_REPEAT);
        nearest_sampler.initialize_handle();
        nearest_sampler.set_param(GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_S, GL_REPEAT);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_T, GL_REPEAT);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_R, GL_REPEAT);

        auto shape_data = wzz::file::read_raw_file_bytes("asset/noiseShapePacked_128_128_128.raw");
        shape_noise.initialize_handle();
        shape_noise.initialize_format_and_data(1, GL_RGBA8, 128, 128, 128,
                                               reinterpret_cast<vec4b *>(shape_data.data()));

        auto detail_data = wzz::file::read_raw_file_bytes("asset/noiseErosionPacked_32_32_32.raw");
        detail_noise.initialize_handle();
        detail_noise.initialize_format_and_data(1, GL_RGBA8, 32, 32, 32, reinterpret_cast<vec4b *>(detail_data.data()));

        auto weather = wzz::image::load_rgba_from_file("asset/weather1.png");
        weather_map.initialize_handle();
        weather_map.initialize_format_and_data(1, GL_RGBA8, weather);

        auto blue = wzz::image::load_rgba_from_file("asset/bluenoise.png");
        blue_noise.initialize_handle();
        blue_noise.initialize_format_and_data(1, GL_RGBA8, blue);

        vao.initialize_handle();

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr, GL_DYNAMIC_DRAW);

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        cloud_buffer.initialize_handle();
        cloud_buffer.reinitialize_buffer_data(nullptr, GL_STATIC_DRAW);

        cloud.fbo.initialize_handle();
        resize(res);
    }

    void resize(const vec2i &res) {
        this->res = res;
        cloud.rbo.destroy();
        cloud.rbo.initialize_handle();
        cloud.rbo.set_format(GL_DEPTH24_STENCIL8, res.x / DownSample, res.y / DownSample);
        cloud.fbo.attach(GL_DEPTH_STENCIL_ATTACHMENT, cloud.rbo);
        cloud.color.destroy();
        cloud.color.initialize_handle();
        cloud.color.initialize_texture(1, GL_RGBA8, res.x / DownSample, res.y / DownSample);
        cloud.fbo.attach(GL_COLOR_ATTACHMENT0, cloud.color);
    }

    void setAtmosphere(const AtmosphereProperties &atmos) {
        atmosphere_buffer.set_buffer_data(&atmos);
    }

    void setCloud(const CloudProperties &cloud) {
        cloud_buffer.set_buffer_data(&cloud);
    }

    void setPerspective(float w_over_h, float camera_fov_rad) {
        params.w_over_h = w_over_h;
        params.scale = tan(camera_fov_rad * 0.5f);
    }

    void setSun(const vec3f &sun_dir) {
        params.sun_dir = sun_dir;
    }

    void render(const texture2d_t &transmittance,
                const texture3d_t &aerial_lut,
                const vec3 &camera_dir,
                const vec3 &camera_right,
                const vec3 &camera_pos) {
        params.view_pos = camera_pos;
        params.camera_dir = camera_dir;
        params.right = camera_right;
        params.up = wzz::math::cross(camera_right, camera_dir).normalized();

        atmosphere_buffer.bind(0);
        cloud_buffer.bind(1);

        params_buffer.set_buffer_data(&params);
        params_buffer.bind(2);

        weather_map.bind(0);
        shape_noise.bind(1);
        detail_noise.bind(2);
        blue_noise.bind(3);
        transmittance.bind(4);
        aerial_lut.bind(5);

        linear_sampler.bind(0);
        linear_sampler.bind(1);
        linear_sampler.bind(2);
        nearest_sampler.bind(3);
        linear_sampler.bind(4);
        linear_sampler.bind(5);

        cloud.fbo.bind();
        cloud.fbo.clear_color_depth_buffer();
        GL_EXPR(glViewport(0, 0, res.x / DownSample, res.y / DownSample));

        shader.bind();
        vao.bind();

        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawArrays(GL_TRIANGLE_STRIP, 0, 4));
        GL_EXPR(glDepthFunc(GL_LESS));

        shader.unbind();
        cloud.fbo.unbind();
        GL_EXPR(glViewport(0, 0, res.x, res.y));
    }

    const texture2d_t &getLUT() const {
        return cloud.color;
    }

private:
    static constexpr int DownSample = 4;
    program_t shader;
    texture3d_t shape_noise;
    texture3d_t detail_noise;
    texture2d_t weather_map;
    texture2d_t blue_noise;
    struct {
        framebuffer_t fbo;
        renderbuffer_t rbo;
        texture2d_t color;
    } cloud;
    vertex_array_t vao;
    vec2i res;

    struct alignas(16) CloudParams {
        vec3 camera_dir;
        float scale;//tan(fov/2)
        vec3 up;
        float w_over_h;
        vec3 right;
        float pad0;
        vec3 view_pos;
        float pad1;
        vec3 sun_dir;
        float pad2;
    } params;
    std140_uniform_block_buffer_t<CloudParams> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    std140_uniform_block_buffer_t<CloudProperties> cloud_buffer;

    sampler_t linear_sampler;
    sampler_t nearest_sampler;
};


class SkyRenderer : public gl_app_t {
public:
    using gl_app_t::gl_app_t;
private:
    void initialize() override {
        // opengl
        GL_EXPR(glEnable(GL_DEPTH_TEST));
        GL_EXPR(glClearColor(0, 0, 0, 0));
        GL_EXPR(glClearDepth(1.0));

        // load model
        loadModel("asset/terrain.obj");

        mesh_shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/mesh.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/mesh.frag"));

        std_unit_atmosphere_properties = preset_atmosphere_properties.toStdUnit();

        transmittance_generator.initialize();
        transmittance_generator.generate(std_unit_atmosphere_properties, transmittance_lut_size);

        multiScattering_generator.initialize();
        multiScattering_generator.setSample(256, 64);
        multiScattering_generator.generate(std_unit_atmosphere_properties, multi_scattering_lut_size,
                                           ground_albedo,
                                           transmittance_generator.getLUT());

        shadow_generator.initialize();

        sky_lut_generator.initialize(sky_lut_size);
        sky_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

        aerial_lut_generator.initialize(aerial_lut_size);
        aerial_lut_generator.setAtmosphere(std_unit_atmosphere_properties);
        aerial_lut_generator.setAerialPerspective(world_scale, max_aerial_distance,
                                                  ray_marching_steps_per_slice,
                                                  enable_shadow, enable_multi_scattering);

        sky_view_renderer.initialize(window->get_window_size());
        sky_view_renderer.enableCloud(enable_cloud);

        mesh_renderer.initialize(window->get_window_size());
        mesh_renderer.setAtmosphere(std_unit_atmosphere_properties);
        mesh_renderer.setAerialPerspective(max_aerial_distance, world_scale);

        sun_disk_renderer.initialize(window->get_window_size());
        sun_disk_renderer.setAtmosphere(std_unit_atmosphere_properties);

        cloud_renderer.initialize(window->get_window_size());
        cloud_renderer.setAtmosphere(std_unit_atmosphere_properties);
        cloud_renderer.setCloud(cloud_properties);
        cloud_renderer.setPerspective(window->get_window_w_over_h(), wzz::math::deg2rad(CameraFovDegree));

        //camera
        camera.set_position({4.087f, 26.7f, 3.957f});
        camera.set_perspective(CameraFovDegree, 0.1f, 100.f);
        camera.set_direction(0, 0.12);
    }

    void frame() override {

        handle_events();

        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 7.f);
        if (ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("Press LCtrl to show/hide cursor");
            ImGui::Text("Use W/A/S/D/Space/LShift to move");
            ImGui::Text("FPS: %.0f", ImGui::GetIO().Framerate);
            if (ImGui::Checkbox("VSync", &vsync)) {
                window->set_vsync(vsync);
            }

            bool update_setting = false;
            update_setting |= ImGui::Checkbox("Enable Multi-Scattering", &enable_multi_scattering);
            update_setting |= ImGui::Checkbox("Enable Shadow", &enable_shadow);
            ImGui::Checkbox("Enable Sky", &enable_sky);
            ImGui::Checkbox("Enable Show Mesh", &enable_draw_mesh);
            ImGui::Checkbox("Enable Sun Disk", &enable_sun_disk);
            if (ImGui::Checkbox("Enable Cloud", &enable_cloud))
                sky_view_renderer.enableCloud(enable_cloud);

            ImGui::Text("Using While Point Color Mapping");
            ImGui::SliderFloat("Exposure", &exposure, 0, 10.f);


            bool update_world_scale = ImGui::InputFloat("World Scale", &world_scale, 1.f);

            ImGui::Text("View Height %.0f in World Scale", camera.get_position().y * world_scale);

            if (ImGui::TreeNode("Atmosphere")) {
                bool update = false;
                update |= ImGui::InputFloat("Planet Radius (km)", &preset_atmosphere_properties.ground_radius, 1.f);
                update |= ImGui::InputFloat("Top Atmosphere Radius (km)",
                                            &preset_atmosphere_properties.top_atmosphere_radius, 1.f);
                update |= ImGui::InputFloat("Ozone Center Altitude (km)", &preset_atmosphere_properties.ozone_center_h,
                                            1.f);
                update |= ImGui::InputFloat("Ozone Thickness (km)", &preset_atmosphere_properties.ozone_width, 1.f);
                update |= ImGui::InputFloat3("Ozone Absorption (10e-6 m^-1)",
                                             &preset_atmosphere_properties.ozone_absorption.x);

                update |= ImGui::InputFloat("Mie Density Height (km)", &preset_atmosphere_properties.mie_density_h,
                                            0.01f);
                update |= ImGui::InputFloat("Mie Absorption (10e-6 m^-1)", &preset_atmosphere_properties.mie_absorption,
                                            0.01f);
                update |= ImGui::InputFloat("Mie Scattering (10e-6 m^-1)", &preset_atmosphere_properties.mie_scattering,
                                            0.01f);
                update |= ImGui::SliderFloat("Mie Asymmetry G", &preset_atmosphere_properties.mie_asymmetry_g, -1.f,
                                             1.f);

                update |= ImGui::InputFloat("Rayleigh Density Height (km)",
                                            &preset_atmosphere_properties.rayleigh_density_h);
                update |= ImGui::InputFloat3("Rayleigh Scattering (10e-6 m^-1)",
                                             &preset_atmosphere_properties.rayleigh_scattering.x);

                update |= ImGui::InputInt2("Transmittance LUT Size", &transmittance_lut_size.x);
                update |= ImGui::InputInt2("Multi-Scattering LUT Size", &multi_scattering_lut_size.x);

                if (update) {
                    std_unit_atmosphere_properties = preset_atmosphere_properties.toStdUnit();

                    transmittance_generator.generate(std_unit_atmosphere_properties, transmittance_lut_size);

                    multiScattering_generator.generate(std_unit_atmosphere_properties, multi_scattering_lut_size,
                                                       ground_albedo,
                                                       transmittance_generator.getLUT());

                    sky_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

                    aerial_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

                    mesh_renderer.setAtmosphere(std_unit_atmosphere_properties);

                    sun_disk_renderer.setAtmosphere(std_unit_atmosphere_properties);

                    cloud_renderer.setAtmosphere(std_unit_atmosphere_properties);
                }
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Sky LUT")) {
                if (ImGui::InputInt2("Sky LUT Size", &sky_lut_size.x)) {
                    sky_lut_generator.resize(sky_lut_size);
                }
                ImGui::InputInt("Sky Ray March Steps", &sky_ray_march_step_count);

                ImGui::TreePop();
            }

            bool update_aerial = false;
            if (ImGui::TreeNode("Aerial LUT")) {
                if (ImGui::InputInt3("Aerial LUT Size", &aerial_lut_size.x)) {
                    aerial_lut_generator.resize(aerial_lut_size);
                }

                update_aerial |= ImGui::SliderFloat("Max Aerial View Distance", &max_aerial_distance, 1000.f, 10000.f);

                update_aerial |= ImGui::InputInt("Ray Marching Steps Per Slice", &ray_marching_steps_per_slice);

                ImGui::TreePop();
            }

            if (update_world_scale || update_aerial || update_setting)
                aerial_lut_generator.setAerialPerspective(world_scale, max_aerial_distance,
                                                          ray_marching_steps_per_slice,
                                                          enable_shadow, enable_multi_scattering);

            if (update_aerial || update_world_scale)
                mesh_renderer.setAerialPerspective(max_aerial_distance, world_scale);

            ImGui::SetNextItemOpen(true, ImGuiCond_Once);
            if (ImGui::TreeNode("Sun")) {
                ImGui::SliderFloat("Sun X Degree", &sun_x_degree, 0, 360);
                ImGui::SliderFloat("Sun Y Degree", &sun_y_degree, -10.f, 80);
                ImGui::InputFloat("Sun Intensity", &sun_intensity);
                ImGui::ColorEdit3("Sun Color", &sun_color.x);
                ImGui::InputFloat("Sun Radius", &sun_radius);
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Volumetric Cloud")) {
                bool update = false;
                update |= ImGui::InputFloat3("Cloud Region Low Corner", &cloud_properties.box_min.x);
                update |= ImGui::InputFloat3("Cloud Region High Corner", &cloud_properties.box_max.x);
                update |= ImGui::InputFloat3("Density to Sigma S", &cloud_properties.density_to_sigma_s.x);
                update |= ImGui::InputFloat3("Density to Sigma T", &cloud_properties.density_to_sigma_t.x);
                update |= ImGui::SliderFloat("Phase G", &cloud_properties.phase_g, -1, 1);
                update |= ImGui::InputInt("Primary Ray March Steps", &cloud_properties.primary_ray_marching_steps);
                update |= ImGui::InputInt("Secondary Ray March Steps", &cloud_properties.secondary_ray_marching_steps);
                update |= ImGui::Checkbox("Enable Multi Scattering",
                                          reinterpret_cast<bool *>(&cloud_properties.enable_multi_scattering));
                update |= ImGui::SliderFloat("Global Coverage", &cloud_properties.g_c, 0, 1);
                update |= ImGui::SliderFloat("Global Density", &cloud_properties.g_d, 0, 1);
                update |= ImGui::SliderFloat("wc0", &cloud_properties.wc0, 0, 1);
                update |= ImGui::SliderFloat("wc1", &cloud_properties.wc1, 0, 1);
                update |= ImGui::SliderFloat("wh", &cloud_properties.wh, 0, 1);
                update |= ImGui::InputFloat("shape tile", &cloud_properties.shape_tile, 0.01);
                update |= ImGui::InputFloat("detail tile", &cloud_properties.detail_tile, 0.01);
                update |= ImGui::SliderFloat("Blend Alpha", &cloud_properties.blend_alpha, 0, 1);
                if (update)
                    cloud_renderer.setCloud(cloud_properties);
                ImGui::TreePop();
            }
        }

        auto [sun_dir, sun_proj_view] = getLight();

        //render
        framebuffer_t::clear_color_depth_buffer();

        shadow_generator.begin();
        for (auto &draw_model: draw_models) {
            shadow_generator.generate(draw_model.vao, draw_model.ebo.index_count(),
                                      draw_model.model, sun_proj_view);
        }
        shadow_generator.end();

        sky_lut_generator.setSun(sun_dir, sun_intensity * sun_color);
        sky_lut_generator.generate(camera.get_position() * world_scale, sky_ray_march_step_count,
                                   enable_multi_scattering,
                                   transmittance_generator.getLUT(),
                                   multiScattering_generator.getLUT());

        aerial_lut_generator.setSun(sun_dir, sun_proj_view);
        aerial_lut_generator.setCamera(camera.get_view_proj().inverse());
        aerial_lut_generator.generate(
                transmittance_generator.getLUT(),
                multiScattering_generator.getLUT(),
                shadow_generator.getShadowMap(),
                camera.get_position());

        framebuffer_t::bind_to_default();
        framebuffer_t::clear_color_depth_buffer();

        if (enable_draw_mesh) {
            mesh_renderer.setSun(sun_dir, sun_intensity * sun_color, sun_proj_view);
            mesh_renderer.begin(transmittance_generator.getLUT(),
                                aerial_lut_generator.getLUT(),
                                shadow_generator.getShadowMap(),
                                camera.get_view_proj());
            for (auto &draw_model: draw_models) {
                mesh_renderer.render(draw_model.model,
                                     draw_model.vao, draw_model.ebo.index_count(), camera.get_position());
            }
            mesh_renderer.end();
        }

        auto camera_dir = camera.get_xyz_direction();
        const vec3f world_up = {0.f, 1.f, 0.f};

        if (enable_cloud) {
            cloud_renderer.setSun(sun_dir);
            cloud_renderer.render(
                    transmittance_generator.getLUT(),
                    aerial_lut_generator.getLUT(),
                    camera_dir,
                    wzz::math::cross(camera_dir, world_up).normalized(),
                    camera.get_position());
        }

        if (enable_sky) {

            sky_view_renderer.render(sky_lut_generator.getLUT(),
                                     cloud_renderer.getLUT(),
                                     camera_dir, wzz::math::cross(camera_dir, world_up).normalized(),
                                     wzz::math::deg2rad(CameraFovDegree), exposure);
        }

        if (enable_sun_disk) {
            auto model = transform::translate(-10.f * sun_dir) * transform::translate(camera.get_position()) *
                         transform::scale(vec3f(sun_radius));
            sun_disk_renderer.setSun(sun_dir, sun_intensity * sun_color);
            sun_disk_renderer.render(transmittance_generator.getLUT(),
                                     multiScattering_generator.getLUT(),
                                     camera.get_view_proj() * model,
                                     camera.get_position().y * world_scale);
        }

        ImGui::End();
        ImGui::PopStyleVar();
    }

    void destroy() override {

    }

private:
    void loadModel(const std::string &filename);

    std::pair<vec3f, mat4> getLight() const;

private:
    AtmosphereProperties preset_atmosphere_properties;
    AtmosphereProperties std_unit_atmosphere_properties;

    CloudProperties cloud_properties;

    static constexpr float CameraFovDegree = 60.f;

    vec2i transmittance_lut_size = {1024, 256};
    vec2i multi_scattering_lut_size = {256, 256};
    vec2i sky_lut_size = {512, 256};
    vec3i aerial_lut_size = {200, 150, 32};

    struct DrawModel {
        vertex_array_t vao;
        vertex_buffer_t<vertex_t> vbo;
        index_buffer_t<uint32_t> ebo;
        mat4 model;
    };
    std::vector<DrawModel> draw_models;

    program_t mesh_shader;

    vec3f ground_albedo = vec3(0.3f);

    struct {
        float sun_x_degree = 0;
        float sun_y_degree = 30;
        float sun_intensity = 1.f;
        vec3f sun_color = vec3(1.f, 1.f, 1.f);
    };

    TransmittanceGenerator transmittance_generator;

    MultiScatteringGenerator multiScattering_generator;

    ShadowGenerator shadow_generator;

    SkyLUTGenerator sky_lut_generator;
    int sky_ray_march_step_count = 40;

    AerialLUTGenerator aerial_lut_generator;
    float max_aerial_distance = 2000.f;
    int ray_marching_steps_per_slice = 4;

    SkyViewRenderer sky_view_renderer;

    MeshRenderer mesh_renderer;

    SunDiskRenderer sun_disk_renderer;
    float sun_radius = 0.07f;

    CloudRenderer cloud_renderer;

    // render resources control
    bool vsync = true;
    bool enable_draw_mesh = true;
    bool enable_multi_scattering = true;
    bool enable_shadow = true;
    bool enable_sky = true;
    bool enable_sun_disk = true;
    bool enable_cloud = true;

    float world_scale = 50.f;
    float exposure = 10.f;

};

void SkyRenderer::loadModel(const std::string &filename) {
    auto model = load_model_from_obj_file(filename);
    for (auto &mesh: model->meshes) {
        draw_models.emplace_back();
        auto &m = draw_models.back();
        m.vao.initialize_handle();
        m.vbo.initialize_handle();
        m.ebo.initialize_handle();
        m.vbo.reinitialize_buffer_data(mesh.vertices.data(), mesh.vertices.size(), GL_STATIC_DRAW);
        m.ebo.reinitialize_buffer_data(mesh.indices.data(), mesh.indices.size(), GL_STATIC_DRAW);
        m.vao.bind_vertex_buffer_to_attrib(attrib_var_t<vec3f>(0), m.vbo, &vertex_t::pos, 0);
        m.vao.bind_vertex_buffer_to_attrib(attrib_var_t<vec3f>(1), m.vbo, &vertex_t::normal, 1);
        m.vao.bind_vertex_buffer_to_attrib(attrib_var_t<vec2f>(2), m.vbo, &vertex_t::tex_coord, 2);
        m.vao.enable_attrib(attrib_var_t<vec3f>(0));
        m.vao.enable_attrib(attrib_var_t<vec3f>(1));
        m.vao.enable_attrib(attrib_var_t<vec2f>(2));
        m.vao.bind_index_buffer(m.ebo);
        m.model = transform::translate(vec3f(0, 15.f, 0));
    }
}

std::pair<vec3f, mat4> SkyRenderer::getLight() const {
    float sun_y_rad = wzz::math::deg2rad(sun_y_degree);
    float sun_x_rad = wzz::math::deg2rad(sun_x_degree);
    float sun_dir_y = std::sin(sun_y_rad);
    float sun_dir_x = std::cos(sun_y_rad) * std::cos(sun_x_rad);
    float sun_dir_z = std::cos(sun_y_rad) * std::sin(sun_x_rad);
    vec3f sun_dir = {sun_dir_x, sun_dir_y, sun_dir_z};
    auto view = transform::look_at(sun_dir * 50.f, {0.f, 0.f, 0.f}, {0.f, 1.f, 0.f});
    auto proj = transform::orthographic(-50.f, 50.f, -50.f, 50.f, 1.f, 200.f);
    return {-sun_dir, proj * view};
}

int main() {
    SkyRenderer(window_desc_t{
            .size = {1600, 900},
            .title = "SkyRenderer",
            .resizeable = false,
            .multisamples = 4,
    }).run();
}