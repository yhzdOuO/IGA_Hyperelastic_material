#ifndef UPDATEDLAGRANGIANELEMENT_H
#define UPDATEDLAGRANGIANELEMENT_H

#include "mechanics/ElementCommon.h"
#include "nurbs/TensorProductNurbsSpace.h"

namespace iga {

/** @brief UL 所需的参考数据；不构造 TL 的 B_A，也不保存参考逆 Jacobian。 */
struct ReferenceGradientData {
    Eigen::VectorXd shapeValues;
    /// 每行是 dN_A/d(xi,eta,zeta)，xi 是 NURBS 参数坐标，不是母单元坐标。
    GradientMatrix parametricGradients;
    /// 每行是 grad_X(N_A)^T，用于计算总变形梯度。
    GradientMatrix referenceGradients;
    /// det(dX/dxi)，只保留标量以供体积比及后续 mean-dilatation 使用。
    double referenceJacobianDeterminant = 0.0;
    /// det(dxi/d(parent))*w_gauss，不含任何物理构形的 Jacobian。
    double parameterIntegrationWeight = 0.0;
};

/** @brief 当前构形梯度、积分测度以及材料点总运动学量。 */
struct SpatialPointData {
    /// 每行是 grad_x(N_A)^T。
    GradientMatrix spatialGradients;
    /// A_t=dx/dxi，参数域到当前构形的映射。
    Eigen::Matrix3d currentJacobian = Eigen::Matrix3d::Zero();
    double currentJacobianDeterminant = 0.0;
    /// J=det(F)，与 det(A_t) 不同。
    double deformationJacobian = 0.0;
    /// 当前体积积分权重 det(A_t)*det(dxi/d(parent))*w_gauss。
    double integrationWeight = 0.0;
    KinematicState kinematics;
};

/** @brief 第一步：根据形函数参数导数和参考控制点计算 grad_X(N)。 */
ReferenceGradientData evaluateReferenceGradients(
    const nurbs::ShapeFunctionData& shape,
    const nurbs::ControlPointMatrix& localReferencePoints,
    double parameterIntegrationWeight);

/**
 * @brief 第二步：根据当前控制点计算 grad_x(N)、dx/dxi、F 和 J。
 * @param localDisplacements 相对初始构形的总位移，不是本次 Newton 增量。
 * @param localCurrentPoints 当前控制点 x_A=X_A+u_A，单元循环外一次构造。
 */
SpatialPointData evaluateSpatialPoint(
    const ReferenceGradientData& reference,
    const nurbs::ControlPointMatrix& localDisplacements,
    const nurbs::ControlPointMatrix& localCurrentPoints);

/**
 * @brief 用当前构形积分单元内力、constitutive 和 initial-stress 刚度。
 *
 * 与 TL 的 computeTotalLagrangianElementResponse() 返回相同的 ElementResult 和自由度布局。
 * 本函数不调用 TL 积分，也不更新初始材料参考构形；目前处理局部材料响应。
 */
ElementResult computeUpdatedLagrangianElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material);

} // namespace iga

#endif // UPDATEDLAGRANGIANELEMENT_H
