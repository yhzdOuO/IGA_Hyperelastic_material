#ifndef ELEMENTCOMMON_H
#define ELEMENTCOMMON_H

#include "mechanics/HyperelasticMaterial.h"
#include "nurbs/ControlPoints.h"

#include <Eigen/Core>

namespace iga {

using GradientMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using BMatrix = Eigen::Matrix<double, 6, 3>;

/** @brief TL/UL 单元向全局装配提供相同的自由度布局和返回类型。 */
struct ElementResult {
    /// (3*n_en) x (3*n_en) 单元一致切线 K_e。
    Eigen::MatrixXd tangentStiffness;
    /// 长度 3*n_en 的单元内力 f_int,e。
    Eigen::VectorXd internalForce;
};

/**
 * @brief 根据一个物理梯度构造 6 x 3 对称梯度算子。
 *
 * 行顺序为 [11,22,33,12,23,13]，剪切行采用工程量。
 * TL 输入 grad_X(N)，UL 输入 grad_x(N)；函数本身与构形无关。
 */
BMatrix makeLinearB(const Eigen::Vector3d& physicalGradient);

/** @brief 由总位移和参考梯度计算 H=du/dX、F=I+H 和 E。 */
KinematicState evaluateKinematics(
    const nurbs::ControlPointMatrix& localDisplacements,
    const GradientMatrix& referenceGradients);

/** @brief 对称应力的 Voigt 向量 [11,22,33,12,23,13]，剪应力不乘 2。 */
Eigen::Matrix<double, 6, 1> stressVoigt(const Eigen::Matrix3d& stress);

} // namespace iga

#endif // ELEMENTCOMMON_H
