#version 330 core
layout(location = 0) in vec3 aPosition;
layout(location = 1) in vec3 aNormal;
layout(location = 2) in float aScalar;

uniform mat4 uMvp;
uniform mat4 uView;
uniform mat3 uNormalMatrix;
out vec3 vNormal;
out vec3 vViewPosition;
out float vScalar;

void main()
{
    // 位置进入裁剪空间，同时保留视空间位置供 headlight 光照使用。
    gl_Position = uMvp * vec4(aPosition, 1.0);
    vViewPosition = (uView * vec4(aPosition, 1.0)).xyz;
    // 法向使用 normal matrix，避免视图矩阵缩放时方向失真。
    vNormal = normalize(uNormalMatrix * aNormal);
    // 标量不在 CPU 转成 RGB，以便运行时调整值域和 colormap。
    vScalar = aScalar;
}
