#version 330 core
uniform vec4 uColor;
out vec4 fragmentColor;
// 当前构型和参考构型通过不同的统一颜色区分。
void main() { fragmentColor = uColor; }
