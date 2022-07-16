#pragma once
#include <CGUtils/api.hpp>
#include <CGUtils/model.hpp>
#include <CGUtils/image.hpp>

#include <cyPoint.h>
#include <cySampleElim.h>


using namespace wzz::gl;
using namespace wzz::model;

constexpr float PI = wzz::math::PI_f;


struct Sphere{
    std::vector<vec3f> positions;
    std::vector<vec3f> normals;
    std::vector<vec2f> uv;
    std::vector<uint32_t> indices;
};

Sphere MakeSphere(int segcount = 64){
    uint32_t X_SEGMENTS = 64;
    uint32_t Y_SEGMENTS = 64;

    Sphere sphere;

    for(uint32_t y = 0; y <= Y_SEGMENTS; y++){
        for(uint32_t x = 0; x <= X_SEGMENTS; x++){
            float x_segment = static_cast<float>(x) / X_SEGMENTS;
            float y_segment = static_cast<float>(y) / Y_SEGMENTS;
            float x_pos = std::cos(x_segment * 2.f * PI) * std::sin(y_segment * PI);
            float y_pos = std::cos(y_segment * PI);
            float z_pos = std::sin(x_segment * 2.f * PI) * std::sin(y_segment * PI);

            sphere.positions.emplace_back(x_pos,y_pos,z_pos);
            sphere.uv.emplace_back(x_segment,y_segment);
            sphere.normals.emplace_back(x_pos,y_pos,z_pos);
        }
    }
    //use GL_TRIANGLE_STRIP
    bool odd_row = false;
    for(uint32_t y = 0; y < Y_SEGMENTS; y++){
        if(!odd_row){
            for(uint32_t x = 0; x <= X_SEGMENTS; x++){
                sphere.indices.emplace_back( y    * (X_SEGMENTS + 1) + x);
                sphere.indices.emplace_back((y+1) * (X_SEGMENTS + 1) + x);
            }
        }
        else{
            for(int x = X_SEGMENTS; x>=0; x--){
                sphere.indices.emplace_back((y + 1) * (X_SEGMENTS + 1) + x);
                sphere.indices.emplace_back( y      * (X_SEGMENTS + 1) + x);
            }
        }
        odd_row = !odd_row;
    }

    return sphere;
}


std::vector<vec2f> getPoissonDiskSamples(int count){
    std::default_random_engine rng{ std::random_device()() };
    std::uniform_real_distribution<float> dis(0, 1);

    std::vector<cy::Point2f> rawPoints;
    for(int i = 0; i < count * 10; ++i)
    {
        const float u = dis(rng);
        const float v = dis(rng);
        rawPoints.push_back({ u, v });
    }

    std::vector<cy::Point2f> outputPoints(count);

    cy::WeightedSampleElimination<cy::Point2f, float, 2> wse;
    wse.SetTiling(true);
    wse.Eliminate(
            rawPoints.data(),    rawPoints.size(),
            outputPoints.data(), outputPoints.size());

    std::vector<vec2f> result;
    for(auto &p : outputPoints)
        result.push_back({ p.x, p.y });

    return result;
}
