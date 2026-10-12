#ifndef TOTALLAGRANGIANELEMENT_H
#define TOTALLAGRANGIANELEMENT_H

#include "mechanics/HyperelasticMaterial.h"
#include "nurbs/TensorProductNurbsSpace.h"

#include <Eigen/Core>
#include <Eigen/LU>
#include <array>
#include <vector>

namespace iga {

using GradientMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3>;
using BMatrix = Eigen::Matrix<double, 6, 3>;

struct ReferencePointData {
    /// n_en x 1 NURBS 基函数值 R_A。
    Eigen::VectorXd shapeValues;
    /// n_en x 3 参数梯度，每行是 dR_A/d(xi,eta,zeta)。
    GradientMatrix parametricGradients;
    /// n_en x 3 参考物理梯度 G_A=grad_X(R_A)。
    GradientMatrix physicalGradients;
    /// 3 x 3 几何 Jacobian dX/dxi。
    Eigen::Matrix3d jacobian = Eigen::Matrix3d::Zero();
    /// 3 x 3 逆几何 Jacobian dxi/dX。
    Eigen::Matrix3d inverseJacobian = Eigen::Matrix3d::Zero();
    /// 几何 Jacobian 行列式，必须为正。
    double jacobianDeterminant = 0.0;
    /// 完整积分权重 det(J)*det(J_parent)*w_gauss。
    double integrationWeight = 0.0;
    /// 所有局部基函数的位移无关 6 x 3 线性 B_A。
    std::vector<BMatrix> linearBMatrices;
};

struct ElementResult {
    /// (3*n_en) x (3*n_en) 单元一致切线 K_e。
    Eigen::MatrixXd tangentStiffness;
    /// 长度 3*n_en 的单元内力 f_int,e。
    Eigen::VectorXd internalForce;
};

/**
 * @brief 构造位移无关的 6 x 3 线性应变-位移矩阵 B_A。
 * @param physicalGradient 参考构形梯度 G_A=grad_X(R_A)。
 * @return 使用工程剪应变顺序 [11,22,33,12,23,13] 的 B_A。
 */
BMatrix makeLinearB(const Eigen::Vector3d& physicalGradient);

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
 * @brief 计算一个积分点上随位移变化的 H、F、E（L2 层）。
 * @param localDisplacements n_en x 3 单元控制变量。
 * @param physicalGradients 关于参考坐标 X 的形函数梯度。
 * @return 当前积分点的纯运动学状态；材料响应由 HyperelasticMaterial 单独计算。
 */
KinematicState evaluateKinematics(
    const nurbs::ControlPointMatrix& localDisplacements,
    const GradientMatrix& physicalGradients);

/**
 * @brief 把对称应力张量转换为固定 Voigt 顺序。
 * @return [S11,S22,S33,S12,S23,S13]^T；剪应力不乘工程因子。
 */
Eigen::Matrix<double, 6, 1> stressVoigt(const Eigen::Matrix3d& stress);

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
ElementResult computeElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material);

} // namespace iga

#endif // TOTALLAGRANGIANELEMENT_H
