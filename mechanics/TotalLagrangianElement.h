#ifndef TOTALLAGRANGIANELEMENT_H
#define TOTALLAGRANGIANELEMENT_H

#include "mechanics/ElementCommon.h"
#include "nurbs/TensorProductNurbsSpace.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <array>
#include <vector>

namespace iga {

struct ReferencePointData {
    /// n_en x 1 NURBS 基函数值 R_A。
    Eigen::VectorXd shapeValues;
    /// n_en x 3 参数梯度，每行是 dR_A/d(xi,eta,zeta)。
    GradientMatrix parametricGradients;
    /// n_en x 3 参考物理梯度 G_A=grad_X(R_A)。
    GradientMatrix referenceGradients;
    /// A_0=dX/dxi，参数域到参考构形的 3 x 3 Jacobian。
    Eigen::Matrix3d referenceJacobian = Eigen::Matrix3d::Zero();
    /// A_0^{-1}=dxi/dX，参考构形到参数域的逆 Jacobian。
    Eigen::Matrix3d inverseReferenceJacobian = Eigen::Matrix3d::Zero();
    /// det(A_0)，必须为正；与变形 Jacobian J=det(F) 不同。
    double referenceJacobianDeterminant = 0.0;
    /// 参考体积积分权重 det(A_0)*det(dxi/d(parent))*w_gauss。
    double integrationWeight = 0.0;
    /// 由 referenceGradients 构造的位移无关 6 x 3 B_A。
    std::vector<BMatrix> referenceBMatrices;
};

/**
 * @brief 计算一个积分点上全部 Total-Lagrangian 参考量（L1 层）。
 * @param gaussIndices 三个方向的 Gauss 点编号 {q_xi,q_eta,q_zeta}。
 * @return 形函数、Jacobian、物理梯度、B_A 和完整积分权重。
 */
ReferencePointData evaluateReferencePoint(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const std::array<int, 3>& gaussIndices);

/**
 * @brief 积分一个单元的一致切线和内力。
 *
 * 本函数刻意止于稠密单元层：不做全局稀疏散射、约束处理或 Newton 迭代，
 * 因而可被静力和动力学装配原样复用。
 *
 * @param space Tensor-product NURBS patch.
 * @param element Element index.
 * @param localDisplacements Element control-point displacements, one point per row.
 * @param material 所有积分点共用的通用超弹性材料模型。
 * @return 稠密单元切线 K_e 和单元内力 f_int,e。
 */
ElementResult computeTotalLagrangianElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material);

} // namespace iga

#endif // TOTALLAGRANGIANELEMENT_H
