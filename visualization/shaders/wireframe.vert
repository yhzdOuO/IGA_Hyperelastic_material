#version 330 core
layout(location = 0) in vec3 aPosition;
uniform mat4 uMvp;
// 线框只需要位置变换，不参与标量着色和光照。
void main() { gl_Position = uMvp * vec4(aPosition, 1.0); }
