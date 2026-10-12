#version 330 core
in vec3 vNormal;
in vec3 vViewPosition;
in float vScalar;
uniform float uMinimum;
uniform float uMaximum;
uniform bool uDiverging;
out vec4 fragmentColor;

vec3 coolwarm(float t)
{
    // 有符号应力：蓝色表示负值、白色表示零附近、红色表示正值。
    if (t < 0.5) return mix(vec3(0.23, 0.30, 0.75), vec3(1.0), 2.0 * t);
    return mix(vec3(1.0), vec3(1.0, 0.20, 0.20), 2.0 * (t - 0.5));
}

vec3 viridis(float t)
{
    // 非负标量的顺序色图近似，保持单调亮度变化。
    return clamp(vec3(0.267 + t * (-0.40 + 1.12 * t),
                      0.005 + t * (1.45 - 0.50 * t),
                      0.329 + t * (0.85 - 1.10 * t)), 0.0, 1.0);
}

void main()
{
    // 用 CPU 两遍统计得到的全局值域将标量归一化到 [0,1]。
    float t = clamp((vScalar - uMinimum) / max(uMaximum - uMinimum, 1e-20), 0.0, 1.0);
    vec3 base = uDiverging ? coolwarm(t) : viridis(t);
    // 相机位置即光源位置；abs 使关闭背面剔除时两面都具有可读光照。
    vec3 lightDirection = normalize(-vViewPosition);
    // 较强环境光保留 colormap 信息，Lambert 漫反射仅增强形状层次。
    float diffuse = abs(dot(normalize(vNormal), lightDirection));
    fragmentColor = vec4(base * (0.42 + 0.58 * diffuse), 1.0);
}
