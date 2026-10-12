#ifndef VISUALIZATION_SURFACEMESH_H
#define VISUALIZATION_SURFACEMESH_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace visualization {

/** @brief 可视化支持的标量场；枚举值同时作为每个顶点标量数组的下标。 */
enum class ScalarField : int {
    AxialStress = 0,
    GravityStress,
    ShearXY,
    ShearYZ,
    ShearXZ,
    VonMises,
    PrincipalMax,
    PrincipalMin,
    Count
};

inline const char* scalarFieldName(ScalarField field)
{
    static const char* names[]{"Axial stress", "Gravity-direction stress",
                               "Shear XY", "Shear YZ", "Shear XZ",
                               "von Mises", "Maximum principal stress",
                               "Minimum principal stress"};
    return names[static_cast<int>(field)];
}

/** @brief 上传到 GPU 的交错顶点布局：位置、解析法向和当前标量值。 */
struct Vertex {
    float position[3]{};
    float normal[3]{};
    float scalar = 0.0f;
};

/**
 * @brief CPU 端的完整表面顶点数据。
 *
 * 参考位置、真实位移及其参数导数始终以 double 保存；draw 只保存当前形变
 * 放大倍数下供 GPU 绘制的 float 数据。这样改变显示倍率不会污染真实应力。
 */
struct SurfaceVertex {
    std::array<double, 3> parameters{};
    Eigen::Vector3d referencePosition = Eigen::Vector3d::Zero();
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero();
    std::array<Eigen::Vector3d, 3> referenceTangents{};
    std::array<Eigen::Vector3d, 3> displacementTangents{};
    std::array<double, static_cast<int>(ScalarField::Count)> scalars{};
    int tangentDirection0 = 0;
    int tangentDirection1 = 1;
    double normalSign = 1.0;
    Vertex draw;
};

/** @brief 边界表面网格及其运行时显示状态。 */
struct SurfaceMesh {
    std::vector<SurfaceVertex> vertices;
    std::vector<std::uint32_t> triangleIndices;
    std::vector<std::uint32_t> wireIndices;
    ScalarField activeField = ScalarField::AxialStress;
    double scalarMinimum = 0.0;
    double scalarMaximum = 0.0;
    double exaggeration = 1.0;

    /** @brief 仅更新绘制位置和法向，不重新计算应力。 */
    void updateGeometry(double scale)
    {
        if (!std::isfinite(scale) || scale <= 0.0) {
            throw std::invalid_argument("deformation scale must be finite and positive");
        }
        exaggeration = scale;
        for (SurfaceVertex& vertex : vertices) {
            // 形变放大只进入 x_draw=X+s*u；应力标量不在这里修改。
            const Eigen::Vector3d position = vertex.referencePosition
                                             + scale * vertex.displacement;
            // 当前构型切向同样采用显示倍率，以保证放大后表面的光照法向正确。
            const Eigen::Vector3d tangent0 =
                vertex.referenceTangents[vertex.tangentDirection0]
                + scale * vertex.displacementTangents[vertex.tangentDirection0];
            const Eigen::Vector3d tangent1 =
                vertex.referenceTangents[vertex.tangentDirection1]
                + scale * vertex.displacementTangents[vertex.tangentDirection1];
            // normalSign 由边界采样器按固定参数面的内/外侧预先确定。
            Eigen::Vector3d normal = vertex.normalSign * tangent0.cross(tangent1);
            const double length = normal.norm();
            if (!std::isfinite(length) || length <= 1.0e-14) {
                throw std::domain_error("degenerate boundary normal");
            }
            normal /= length;
            for (int component = 0; component < 3; ++component) {
                vertex.draw.position[component] = static_cast<float>(position(component));
                vertex.draw.normal[component] = static_cast<float>(normal(component));
            }
        }
    }

    /** @brief 切换标量字段，并完成上传前的全局值域统计。 */
    void selectField(ScalarField field)
    {
        activeField = field;
        scalarMinimum = std::numeric_limits<double>::infinity();
        scalarMaximum = -std::numeric_limits<double>::infinity();
        const int index = static_cast<int>(field);
        // 第一遍选择标量并统计实际的全局极值。
        for (SurfaceVertex& vertex : vertices) {
            const double value = vertex.scalars[static_cast<std::size_t>(index)];
            if (!std::isfinite(value)) {
                throw std::domain_error("non-finite visualization scalar");
            }
            vertex.draw.scalar = static_cast<float>(value);
            scalarMinimum = std::min(scalarMinimum, value);
            scalarMaximum = std::max(scalarMaximum, value);
        }
        if (vertices.empty()) {
            scalarMinimum = scalarMaximum = 0.0;
        }
        // 有符号应力采用关于零对称的值域，使拉压两侧颜色可直接比较。
        if (field != ScalarField::VonMises) {
            const double magnitude = std::max(std::abs(scalarMinimum),
                                              std::abs(scalarMaximum));
            scalarMinimum = -magnitude;
            scalarMaximum = magnitude;
        } else {
            scalarMinimum = 0.0;
        }
        // 常量场也必须给 shader 一个非零区间，避免归一化除零。
        if (scalarMaximum - scalarMinimum < 1.0e-14) {
            scalarMaximum = scalarMinimum + 1.0;
        }
    }

    /** @brief von Mises 非负，使用顺序色图；其余字段默认使用发散色图。 */
    bool usesDivergingColorMap() const { return activeField != ScalarField::VonMises; }
};

} // namespace visualization

#endif // VISUALIZATION_SURFACEMESH_H
