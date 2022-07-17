#version 460 core

layout(location = 0) in vec3 iVertexPos;
layout(location = 1) in vec3 iVertexNormal;
layout(location = 2) in vec2 iTexCoord;

uniform mat4 Model;
uniform mat4 ProjView;

layout(location = 0) out vec3 oVertexPos;
layout(location = 1) out vec3 oVertexNormal;
layout(location = 2) out vec2 oTexCoord;
layout(location = 3) out vec4 oScreenCoord;

void main() {
    gl_Position = ProjView * Model * (vec4(iVertexPos, 1.0));
    oVertexPos = vec3(Model * vec4(iVertexPos, 1.f));
    //ok for simple scale or rotate or transform
    oVertexNormal = vec3(Model * vec4(iVertexNormal, 0.0));
    oTexCoord = iTexCoord;

    oScreenCoord = gl_Position;// -1 ~ 1 for xyz/w
}
