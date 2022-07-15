#include "atmosphere.hpp"

inline vec3i getGroupSize(int x,int y = 1,int z = 1){
    constexpr int group_thread_size_x = 16;
    constexpr int group_thread_size_y = 16;
    constexpr int group_thread_size_z = 16;
    const int group_size_x = (x + group_thread_size_x - 1) / group_thread_size_x;
    const int group_size_y = (y + group_thread_size_y - 1) / group_thread_size_y;
    const int group_size_z = (z + group_thread_size_z - 1) / group_thread_size_z;
    return {group_size_x,group_size_y,group_size_z};
}

class TransmittanceGenerator{
public:
    void initialize(){
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/transmittance.comp"));
        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr,GL_STATIC_DRAW);
    }

    void generate(const AtmosphereProperties& ap,const vec2i& lut_size){
        atmosphere_buffer.set_buffer_data(&ap);
        atmosphere_buffer.bind(0);

        if(lut.handle()) lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1,GL_RGBA32F,lut_size.x,lut_size.y);

        lut.bind_image(0,0,GL_READ_WRITE,GL_RGBA32F);
        shader.bind();

        auto group_size = getGroupSize(lut_size.x,lut_size.y);
        GL_EXPR(glDispatchCompute(group_size.x,group_size.y,1));

        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));

        shader.unbind();
    }
    const texture2d_t& getLUT() const{
        return lut;
    }
private:
    program_t shader;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    texture2d_t lut;
};

class MultiScatteringGenerator{
public:
    void initialize(){
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/multiscattering.comp"));
        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr,GL_STATIC_DRAW);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(&params,GL_STATIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
    }
    void setSample(int ray_march_steps,int dir_samples){
        params.ray_march_step_count = ray_march_steps;
        if(params.dir_sample_count != dir_samples){
            params.dir_sample_count = dir_samples;
            dir_samples_buffer.destroy();
            dir_samples_buffer.initialize_handle();
            auto dir_samples_arr = getPoissonDiskSamples(dir_samples);
            dir_samples_buffer.initialize_buffer_data(dir_samples_arr.data(),dir_samples_arr.size() * sizeof(vec2f),GL_DYNAMIC_STORAGE_BIT);
        }
    }
    void generate(const AtmosphereProperties& ap,const vec2i& lut_size,
                  const vec3f& albedo,
                  const texture2d_t& transmittance){
        atmosphere_buffer.set_buffer_data(&ap);
        atmosphere_buffer.bind(0);

        lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1,GL_RGBA32F,lut_size.x,lut_size.y);
        lut.bind_image(0,0,GL_READ_WRITE,GL_RGBA32F);

        params.ground_albedo = albedo;
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(1);

        dir_samples_buffer.bind(0);

        transmittance.bind(0);
        linear_sampler.bind(0);

        shader.bind();

        auto group_size = getGroupSize(lut_size.x,lut_size.y);
        GL_EXPR(glDispatchCompute(group_size.x,group_size.y,1));

        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));

        shader.unbind();
    }
    const texture2d_t& getLUT() const{
        return lut;
    }
private:
    struct Params{
        vec3f ground_albedo = vec3(0.3f);
        int dir_sample_count = 0;
        vec3f sun_intensity = vec3f(1);//todo
        int ray_march_step_count = 0;
    }params;
    std140_uniform_block_buffer_t<Params> params_buffer;
    storage_buffer_t<vec2f> dir_samples_buffer;
    program_t shader;
    texture2d_t lut;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    sampler_t linear_sampler;
};

class ShadowGenerator{
public:
    void initialize(){
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/mesh.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/depth.frag"));
        shadow.fbo.initialize_handle();
        shadow.rbo.initialize_handle();
        shadow.rbo.set_format(GL_DEPTH32F_STENCIL8,4096,4096);
        shadow.fbo.attach(GL_DEPTH_STENCIL_ATTACHMENT,shadow.rbo);
        shadow.tex.initialize_handle();
        shadow.tex.initialize_texture(1,GL_R32F,4096,4096);
        shadow.fbo.attach(GL_COLOR_ATTACHMENT0,shadow.tex);
        assert(shadow.fbo.is_complete());
    }
    void begin(){
        shadow.fbo.bind();
        GL_EXPR(glViewport(0,0,4096,4096));
        GL_EXPR(glDrawBuffer(GL_COLOR_ATTACHMENT0));
        shadow.fbo.clear_color_depth_buffer();
        shader.bind();
    }
    void generate(const vertex_array_t& vao,GLsizei count,const mat4& model,const mat4& proj_view){

        shader.set_uniform_var("Model",model);
        shader.set_uniform_var("ProjView",proj_view);
        vao.bind();
        GL_EXPR(glDrawElements(GL_TRIANGLES,count,GL_UNSIGNED_INT,nullptr));
        vao.unbind();

    }
    void end(){
        shader.unbind();
        shadow.fbo.unbind();
    }
    const texture2d_t& getShadowMap() const{
        return shadow.tex;
    }
private:
    program_t shader;
    struct {
        framebuffer_t fbo;
        renderbuffer_t rbo;
        texture2d_t tex;
    }shadow;
};

class SkyLUTGenerator{
public:
    void initialize(const vec2i& res){
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/quad.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/sky_lut.frag"));

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr,GL_STATIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);

        sky.fbo.initialize_handle();
        resize(res);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr,GL_DYNAMIC_DRAW);

        vao.initialize_handle();
    }
    void resize(const vec2i& res){
        sky.rbo.destroy();
        sky.rbo.initialize_handle();
        sky.rbo.set_format(GL_DEPTH24_STENCIL8,res.x,res.y);
        sky.fbo.attach(GL_DEPTH_STENCIL_ATTACHMENT,sky.rbo);
        sky.lut.destroy();
        sky.lut.initialize_handle();
        sky.lut.initialize_texture(1,GL_RGBA32F,res.x,res.y);
        sky.fbo.attach(GL_COLOR_ATTACHMENT0,sky.lut);
        assert(sky.fbo.is_complete());
        sky_res = res;
    }
    void setAtmosphere(const AtmosphereProperties& atmos){
        atmosphere_buffer.set_buffer_data(&atmos);
    }
    void setSun(const vec3f& sun_dir,const vec3f& sun_radiance){
        params.sun_dir = sun_dir;
        params.sun_radiance = sun_radiance;
    }
    void generate(const vec3f& view_pos,int ray_march_steps,bool enable_multi_scattering,
                  const texture2d_t& transmittance,
                  const texture2d_t& multi_scattering){
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
        GL_EXPR(glViewport(0,0,sky_res.x,sky_res.y));
        GL_EXPR(glDrawBuffer(GL_COLOR_ATTACHMENT0));
        sky.fbo.clear_color_depth_buffer();
        shader.bind();
        vao.bind();

        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawArrays(GL_TRIANGLE_STRIP,0,4));
        GL_EXPR(glDepthFunc(GL_LESS));

        vao.unbind();
        shader.unbind();
        sky.fbo.unbind();
    }
    const texture2d_t& getLUT() const{
        return sky.lut;
    }
private:
    struct Params{
        vec3f view_pos;
        int ray_march_step_count;
        vec3f sun_dir;
        int enable_multi_scattering;
        vec3 sun_radiance;
        float pad = 0.f;
    }params;
    std140_uniform_block_buffer_t<Params> params_buffer;
    program_t shader;
    sampler_t linear_sampler;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    struct{
        framebuffer_t fbo;
        renderbuffer_t rbo;
        texture2d_t lut;
    }sky;
    vertex_array_t vao;
    vec2i sky_res;
};

class AerialLUTGenerator{
public:
    void initialize(const vec3i& res){
        shader = program_t::build_from(
                shader_t<GL_COMPUTE_SHADER>::from_file("asset/glsl/aerial_lut.comp"));

        atmosphere_buffer.initialize_handle();
        atmosphere_buffer.reinitialize_buffer_data(nullptr,GL_STATIC_DRAW);

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr,GL_DYNAMIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);
        nearest_sampler.initialize_handle();
        nearest_sampler.set_param(GL_TEXTURE_MIN_FILTER,GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_MAG_FILTER,GL_NEAREST);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        nearest_sampler.set_param(GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);

        resize(res);

    }
    void resize(const vec3i& res){
        lut.destroy();
        lut.initialize_handle();
        lut.initialize_texture(1,GL_RGBA32F,res.x,res.y,res.z);
        lut_size = res;

        params.slice_z_count = res.z;

    }
    void setAtmosphere(const AtmosphereProperties& atmos){
        atmosphere_buffer.set_buffer_data(&atmos);
    }
    void setCamera(const mat4& inv_proj_view){
        const vec4f A0 = inv_proj_view * vec4f(-1, 1, 0.2f, 1);
        const vec4f A1 = inv_proj_view *vec4f(-1, 1, 0.5f, 1);

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
    void setSun(const vec3& sun_dir,const mat4& sun_proj_view){
        params.sun_dir = sun_dir;
        params.sun_theta = std::asin(-sun_dir.y);
        params.shadow_vp = sun_proj_view;
    }
    void setAerialPerspective(float world_scale,
                              float max_aerial_distance,
                              int ray_march_steps_per_slice,
                              bool enable_shadow,
                              bool enable_multi_scattering){
        params.world_scale = world_scale;
        params.max_view_distance = max_aerial_distance;
        params.ray_march_step_count_per_slice = ray_march_steps_per_slice;
        params.enable_shadow = static_cast<int>(enable_shadow);
        params.enable_multiscattering = static_cast<int>(enable_multi_scattering);
    }
    void generate(
            const texture2d_t& transmittance,
            const texture2d_t& multi_scattering,
            const texture2d_t& shadow_map,
            const vec3f& view_pos
            ){

        atmosphere_buffer.bind(0);

        params.view_pos = view_pos;
        params.view_height = view_pos.y * params.world_scale;
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(1);

        transmittance.bind(0);
        transmittance.bind(1);
        shadow_map.bind(2);
        linear_sampler.bind(0);
        linear_sampler.bind(1);
        nearest_sampler.bind(2);

        lut.bind_image(0,0,GL_READ_WRITE,GL_RGBA32F);

        shader.bind();
        auto group_size = getGroupSize(lut_size.x,lut_size.y,1);
        GL_EXPR(glDispatchCompute(group_size.x,group_size.y,1));
        GL_EXPR(glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT));
        shader.unbind();

    }
private:
    struct AerialParams{
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
    }params;
    std140_uniform_block_buffer_t<AerialParams> params_buffer;
    std140_uniform_block_buffer_t<AtmosphereProperties> atmosphere_buffer;
    texture3d_t lut;
    vec3i lut_size;
    sampler_t linear_sampler;
    sampler_t nearest_sampler;
    program_t shader;
};

class MeshRenderer{
public:

private:

};

class SkyViewRenderer{
public:
    void initialize(const vec2i& res){
        shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/quad.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/sky.frag"));

        vao.initialize_handle();

        params_buffer.initialize_handle();
        params_buffer.reinitialize_buffer_data(nullptr,GL_DYNAMIC_DRAW);

        linear_sampler.initialize_handle();
        linear_sampler.set_param(GL_TEXTURE_MIN_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_MAG_FILTER,GL_LINEAR);
        linear_sampler.set_param(GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
        linear_sampler.set_param(GL_TEXTURE_WRAP_R,GL_CLAMP_TO_EDGE);

        resize(res);
    }
    void resize(const vec2i& res){
        view_res = res;
        shader.bind();
        shader.set_uniform_var("WOverH",res.x * 1.f / res.y);
        shader.unbind();
    }
    void render(const texture2d_t& sky_lut,
                const vec3& camera_dir,
                const vec3& camera_right,
                float camera_fov_rad,
                float exposure,
                bool enable_tone){
        params.camera_dir = camera_dir;
        params.right = camera_right;
        params.up = wzz::math::cross(camera_right,camera_dir).normalized();
        params.exposure = exposure;
        params.enable_tone_mapping = static_cast<int>(enable_tone);
        params.scale = std::tan(0.5f * camera_fov_rad);
        params_buffer.set_buffer_data(&params);
        params_buffer.bind(0);

        framebuffer_t::clear_color_depth_buffer();
        GL_EXPR(glViewport(0,0,view_res.x,view_res.y));
        shader.bind();
        sky_lut.bind(0);
        linear_sampler.bind(0);
        vao.bind();

        GL_EXPR(glDepthFunc(GL_LEQUAL));
        GL_EXPR(glDrawArrays(GL_TRIANGLE_STRIP,0,4));
        GL_EXPR(glDepthFunc(GL_LESS));

        vao.unbind();
        shader.unbind();

    }
private:
    struct SkyParams{
        vec3f camera_dir;
        float exposure;
        vec3 up;
        int enable_tone_mapping;
        vec3 right;
        float scale;//tan(fov/2)
    }params;
    std140_uniform_block_buffer_t<SkyParams> params_buffer;
    program_t shader;
    vertex_array_t vao;
    vec2i view_res;
    sampler_t linear_sampler;
};

class SunDiskRenderer{
public:

private:

};


class SkyRenderer:public gl_app_t{
public:
    using gl_app_t::gl_app_t;
private:
    void initialize() override{
        // opengl
        GL_EXPR(glEnable(GL_DEPTH_TEST));
        GL_EXPR(glClearColor(0, 0, 0, 0));
        GL_EXPR(glClearDepth(1.0));



        // load model
        loadModel("asset/terrain5.obj");


        mesh_shader = program_t::build_from(
                shader_t<GL_VERTEX_SHADER>::from_file("asset/glsl/mesh.vert"),
                shader_t<GL_FRAGMENT_SHADER>::from_file("asset/glsl/mesh.frag"));

        std_unit_atmosphere_properties = preset_atmosphere_properties.toStdUnit();

        transmittance_generator.initialize();
        transmittance_generator.generate(std_unit_atmosphere_properties,transmittance_lut_size);

        multiScattering_generator.initialize();
        multiScattering_generator.setSample(256,64);
        multiScattering_generator.generate(std_unit_atmosphere_properties,multi_scattering_lut_size,
                                           ground_albedo,
                                           transmittance_generator.getLUT());

        shadow_generator.initialize();

        sky_lut_generator.initialize(sky_lut_size);
        sky_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

        aerial_lut_generator.initialize(aerial_lut_size);
        aerial_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

        sky_view_renderer.initialize(window->get_window_size());

        //camera
        camera.set_position({4.087f,3.7f,3.957f});
        camera.set_perspective(CameraFovDegree,0.1f,100.f);
        camera.set_direction(-PI,0);
    }

    void frame() override{

        handle_events();

        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding,7.f);
        if (ImGui::Begin("Settings", nullptr,ImGuiWindowFlags_AlwaysAutoResize)){
            ImGui::Text("Press LCtrl to show/hide cursor");
            ImGui::Text("Use W/A/S/D/Space/LShift to move");
            ImGui::Text("FPS: %.0f", ImGui::GetIO().Framerate);
            if(ImGui::Checkbox("VSync",&vsync)){
                window->set_vsync(vsync);
            }


            ImGui::Checkbox("Enable Multi-Scattering",&enable_multi_scattering);
            ImGui::Checkbox("Enable Shadow",&enable_shadow);
            ImGui::Checkbox("Enable Sky",&enable_sky);
            ImGui::Checkbox("Enable Show Mesh",&enable_draw_mesh);
            ImGui::Checkbox("Enable Sun Disk",&enable_sun_disk);
            ImGui::Checkbox("Enable Tone Mapping",&enable_tone_mapping);
            if(!enable_tone_mapping){
                ImGui::Text("Using While Point Color Mapping");
                ImGui::SliderFloat("Exposure",&exposure,0,10.f);
            }

            bool update_world_scale = ImGui::InputFloat("World Scale",&world_scale,1.f);

            if(ImGui::TreeNode("Atmosphere")){
                bool update = false;
                update |= ImGui::InputFloat("Planet Radius (km)",&preset_atmosphere_properties.ground_radius,1.f);
                update |= ImGui::InputFloat("Top Atmosphere Radius (km)",&preset_atmosphere_properties.top_atmosphere_radius,1.f);
                update |= ImGui::InputFloat("Ozone Center Altitude (km)",&preset_atmosphere_properties.ozone_center_h,1.f);
                update |= ImGui::InputFloat("Ozone Thickness (km)",&preset_atmosphere_properties.ozone_width,1.f);
                update |= ImGui::InputFloat3("Ozone Absorption (10e-6 m^-1)",&preset_atmosphere_properties.ozone_absorption.x);

                update |= ImGui::InputFloat("Mie Density Height (km)",&preset_atmosphere_properties.mie_density_h,0.01f);
                update |= ImGui::InputFloat("Mie Absorption (10e-6 m^-1)",&preset_atmosphere_properties.mie_absorption,0.01f);
                update |= ImGui::InputFloat("Mie Scattering (10e-6 m^-1)",&preset_atmosphere_properties.mie_scattering,0.01f);
                update |= ImGui::SliderFloat("Mie Asymmetry G",&preset_atmosphere_properties.mie_asymmetry_g,-1.f,1.f);

                update |= ImGui::InputFloat("Rayleigh Density Height (km)",&preset_atmosphere_properties.rayleigh_density_h);
                update |= ImGui::InputFloat3("Rayleigh Scattering (10e-6 m^-1)",&preset_atmosphere_properties.rayleigh_scattering.x);

                update |= ImGui::InputInt2("Transmittance LUT Size",&transmittance_lut_size.x);
                update |= ImGui::InputInt2("Multi-Scattering LUT Size",&multi_scattering_lut_size.x);

                if(update){
                    std_unit_atmosphere_properties = preset_atmosphere_properties.toStdUnit();
                    transmittance_generator.generate(std_unit_atmosphere_properties,transmittance_lut_size);

                    multiScattering_generator.generate(std_unit_atmosphere_properties,multi_scattering_lut_size,
                                                       ground_albedo,
                                                       transmittance_generator.getLUT());

                    sky_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

                    aerial_lut_generator.setAtmosphere(std_unit_atmosphere_properties);

                }
                ImGui::TreePop();
            }
            if(ImGui::TreeNode("Sky LUT")){
                if(ImGui::InputInt2("Sky LUT Size",&sky_lut_size.x)){
                    sky_lut_generator.resize(sky_lut_size);
                }
                ImGui::InputInt("Sky Ray March Steps",&sky_ray_march_step_count);

                ImGui::TreePop();
            }

            bool update_aerial = false;
            if(ImGui::TreeNode("Aerial LUT")){
                if(ImGui::InputInt3("Aerial LUT Size",&aerial_lut_size.x)){
                    aerial_lut_generator.resize(aerial_lut_size);
                }

                update_aerial |= ImGui::SliderFloat("Max Aerial View Distance",&max_aerial_distance,1000.f,10000.f);

                update_aerial |= ImGui::InputInt("Ray Marching Steps Per Slice",&ray_marching_steps_per_slice);

                ImGui::TreePop();
            }
            if(update_world_scale || update_aerial)
                aerial_lut_generator.setAerialPerspective(world_scale,max_aerial_distance,ray_marching_steps_per_slice,
                                                      enable_shadow,enable_multi_scattering);

            if(ImGui::TreeNode("Sun")){
                ImGui::SliderFloat("Sun X Degree",&sun_x_degree,0,360);
                ImGui::SliderFloat("Sun Y Degree",&sun_y_degree,0,90);
                ImGui::InputFloat3("Sun Radiance",&sun_radiance.x);
                ImGui::TreePop();
            }

            if(ImGui::TreeNode("Volumetric Cloud")){

                ImGui::TreePop();
            }
        }

        auto [sun_dir,sun_proj_view] = getLight();


        //render
        framebuffer_t::clear_color_depth_buffer();


        shadow_generator.begin();
        for(auto& draw_model:draw_models){
            shadow_generator.generate(draw_model.vao,draw_model.ebo.index_count(),
                                      draw_model.model,sun_proj_view);
        }
        shadow_generator.end();


        sky_lut_generator.setSun(sun_dir,sun_radiance);
        sky_lut_generator.generate(camera.get_position() * world_scale
                                   ,sky_ray_march_step_count,
                                   enable_multi_scattering,
                                   transmittance_generator.getLUT(),
                                   multiScattering_generator.getLUT());

        aerial_lut_generator.setSun(sun_dir,sun_proj_view);
        aerial_lut_generator.setCamera(camera.get_view_proj().inverse());
        aerial_lut_generator.generate(
                transmittance_generator.getLUT(),
                multiScattering_generator.getLUT(),
                shadow_generator.getShadowMap(),
                camera.get_position());

        auto camera_dir = camera.get_xyz_direction();
        const vec3f world_up = {0.f,1.f,0.f};

        if(enable_sky){
            sky_view_renderer.render(sky_lut_generator.getLUT(),
                                     camera_dir,wzz::math::cross(camera_dir,world_up).normalized(),
                                     wzz::math::deg2rad(CameraFovDegree),exposure,enable_tone_mapping);
        }

        if(enable_draw_mesh){
            mesh_shader.bind();
            mesh_shader.set_uniform_var("ProjView",camera.get_view_proj());
            for(auto& draw_model:draw_models){
                draw_model.vao.bind();
                mesh_shader.set_uniform_var("Model",draw_model.model);

                GL_EXPR(glDrawElements(GL_TRIANGLES,draw_model.ebo.index_count(),GL_UNSIGNED_INT,nullptr));

                draw_model.vao.unbind();
            }
            mesh_shader.unbind();
        }

        if(enable_sky){

        }


        if(enable_sun_disk){

        }

        ImGui::End();
        ImGui::PopStyleVar();
    }

    void destroy() override{

    }

private:
    void loadModel(const std::string& filename);

    std::pair<vec3f,mat4> getLight() const;
private:
    AtmosphereProperties preset_atmosphere_properties;
    AtmosphereProperties std_unit_atmosphere_properties;

    static constexpr float CameraFovDegree = 60.f;

    vec2i transmittance_lut_size = {256,256};
    vec2i multi_scattering_lut_size = {256,256};
    vec2i sky_lut_size = {192,108};
    vec3i aerial_lut_size = {200,150,32};


    struct DrawModel{
        vertex_array_t vao;
        vertex_buffer_t<vertex_t> vbo;
        index_buffer_t<uint32_t> ebo;
        mat4 model;
    };
    std::vector<DrawModel> draw_models;

    program_t mesh_shader;

    vec3f ground_albedo = vec3(0.3f);

    struct{
        float sun_x_degree = 0;
        float sun_y_degree = 60;
        vec3f sun_radiance = {10.f,10.f,10.f};
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

    // render resources control
    bool vsync = true;
    bool enable_draw_mesh = true;
    bool enable_multi_scattering = true;
    bool enable_shadow = true;
    bool enable_sky = true;
    bool enable_sun_disk = true;
    bool enable_tone_mapping = false;

    float world_scale = 200.f;
    float exposure = 2.f;

};

void SkyRenderer::loadModel(const std::string &filename) {
    auto model = load_model_from_obj_file(filename);
    for(auto& mesh:model->meshes){
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

    }
}

std::pair<vec3f, mat4> SkyRenderer::getLight() const{
    float sun_y_rad = wzz::math::deg2rad(sun_y_degree);
    float sun_x_rad = wzz::math::deg2rad(sun_x_degree);
    float sun_dir_y = std::sin(sun_y_rad);
    float sun_dir_x = std::cos(sun_y_rad) * std::cos(sun_x_rad);
    float sun_dir_z = std::cos(sun_y_rad) * std::sin(sun_x_rad);
    vec3f sun_dir = {sun_dir_x,sun_dir_y,sun_dir_z};
    auto view = transform::look_at(sun_dir * 50.f,{0.f,0.f,0.f},{0.f,1.f,0.f});
    auto proj = transform::orthographic(-50.f,50.f,-50.f,50.f,1.f,200.f);
    return {-sun_dir,proj * view};
}

int main(){
    SkyRenderer(window_desc_t{
        .size = {1280,720},
        .title = "SkyRenderer",
        .resizeable = false,
        .multisamples = 4,
    }).run();
}