#ifndef VISUALIZATION_FIELDEVALUATOR_H
#define VISUALIZATION_FIELDEVALUATOR_H

#include "mechanics/TotalLagrangianElement.h"
#include "visualization/SurfaceMesh.h"

#include <Eigen/Core>
#include <array>

namespace visualization {

/** @brief 一个任意参数点上的几何、位移、切向和应力求值结果。 */
struct FieldPoint {
    Eigen::Vector3d referencePosition = Eigen::Vector3d::Zero();
    Eigen::Vector3d displacement = Eigen::Vector3d::Zero();
    std::array<Eigen::Vector3d, 3> referenceTangents{};
    std::array<Eigen::Vector3d, 3> displacementTangents{};
    Eigen::Matrix3d cauchyStress = Eigen::Matrix3d::Zero();
    std::array<double, static_cast<int>(ScalarField::Count)> scalars{};
};

/**
 * @brief 在任意参数坐标处复用求解器公式计算位移和 Cauchy 应力。
 *
 * 本类不负责采样和渲染，只负责把 patch、控制点位移与材料响应转换成
 * 单点物理量，因此可以独立进行数值测试。
 */
class FieldEvaluator {
public:
    /** @brief 保存求值所需输入，并将轴向、重力方向归一化。 */
    FieldEvaluator(const nurbs::TensorProductNurbsSpace& space,
                   const Eigen::VectorXd& globalDisplacements,
                   const iga::HyperelasticMaterial& material,
                   const Eigen::Vector3d& axialDirection,
                   const Eigen::Vector3d& gravityDirection);

    /** @brief 在 (xi,eta,zeta) 处执行 NURBS、运动学和应力求值。 */
    FieldPoint evaluate(const std::array<double, 3>& parameters) const;

private:
    const nurbs::TensorProductNurbsSpace& space_;
    Eigen::VectorXd displacements_;
    /// 材料由调用方持有，并且必须覆盖本求值器的使用期。
    const iga::HyperelasticMaterial& material_;
    Eigen::Vector3d axialDirection_;
    Eigen::Vector3d gravityDirection_;
};

} // namespace visualization

#endif // VISUALIZATION_FIELDEVALUATOR_H
