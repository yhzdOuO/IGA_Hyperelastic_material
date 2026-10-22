#include "mechanics/TotalLagrangianElement.h"

#include <cmath>
#include <stdexcept>

namespace iga {

ReferencePointData evaluateReferencePoint(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const std::array<int, 3>& gaussIndices)
{
    // ---------- L1：参考构形中的几何量（Newton 迭代中保持不变） ----------
    const nurbs::ShapeFunctionData shape =
        space.evaluateShapeFunctions(element, gaussIndices); // &计算型函数和参数梯度: R 和 \nabla R
    const nurbs::ControlPointMatrix localReferencePoints =
        space.elementControlPoints(element);

    ReferencePointData result;
    result.shapeValues = shape.values;
    result.parametricGradients = shape.parametricGradients;
    // X_e: n_en x 3，dN/dxi: n_en x 3，因此
    // A_0 = X_e^T * dN/dxi 的大小为 3x3，(A_0)_ij=dX_i/dxi_j。
    result.referenceJacobian.noalias() = localReferencePoints.transpose()
                                         * result.parametricGradients;
    result.referenceJacobianDeterminant = result.referenceJacobian.determinant();
    if (!std::isfinite(result.referenceJacobianDeterminant)
        || result.referenceJacobianDeterminant <= 0.0) {
        throw std::domain_error(
            "evaluateReferencePoint: geometry Jacobian determinant must be positive");
    }
    result.inverseReferenceJacobian = result.referenceJacobian.inverse();
    // dN/dX = dN/dxi * A_0^{-1}，结果仍为 n_en x 3；梯度按行保存。
    result.referenceGradients.noalias() = result.parametricGradients
                                          * result.inverseReferenceJacobian;
    // 积分权重同时包含几何 Jacobian、母单元 Jacobian 和三维 Gauss 权重。
    result.integrationWeight = result.referenceJacobianDeterminant
        * space.parentJacobianDeterminant(element)
        * space.gaussWeight(gaussIndices); //& total lagrangian下，当前构型改变不用更新Jacobian，该Jacobian是参数空间到参考构型的映射，并非到当前构型

    // 每个局部基函数保存一个 6x3 的 B_A，后续各 Newton 迭代可复用。
    result.referenceBMatrices.reserve(static_cast<std::size_t>(space.numElementNodes()));
    for (int a = 0; a < space.numElementNodes(); ++a) {
        const Eigen::Vector3d gradient = result.referenceGradients.row(a).transpose();
        result.referenceBMatrices.push_back(makeLinearB(gradient));
    }
    return result; //型函数+参数梯度，参数空间到物理空间的Jacobian矩阵及其行列式，物理梯度，积分权重
}

ElementResult computeTotalLagrangianElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material)
{
    const int nodes = space.numElementNodes();
    if (localDisplacements.rows() != nodes || localDisplacements.cols() != 3
        || !localDisplacements.allFinite()) {
        throw std::invalid_argument(
            "computeTotalLagrangianElementResponse: local displacement dimensions are inconsistent");
    }

    // Ke 为 (3*n_en) x (3*n_en)，fe 为 (3*n_en) x 1；
    // 局部自由度按 [u0x,u0y,u0z,u1x,...] 排列。
    ElementResult result;
    result.tangentStiffness = Eigen::MatrixXd::Zero(3 * nodes, 3 * nodes);
    result.internalForce = Eigen::VectorXd::Zero(3 * nodes);
    // ---------- 三维 Gauss 积分循环 ----------
    for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
        for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
            for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                const std::array<int, 3> gaussIndices{q0, q1, q2};
                // 1. N、参数梯度和参考梯度，以及参考积分权重。
                const ReferencePointData reference =
                    evaluateReferencePoint(space, element, gaussIndices); // 计算型函数，NURBS相关
                // 2. 总运动学量 H、F、E。
                const KinematicState kinematics = evaluateKinematics(
                    localDisplacements, reference.referenceGradients); // 位移d相关的运动学量
                // 3. 材料提供第二 PK 应力和材料弹性矩阵。
                const MaterialResponse response = material.evaluate(kinematics);
                if (!response.secondPiolaStress.allFinite()
                    || !response.materialTangent.allFinite()
                    || !std::isfinite(response.strainEnergyDensity)) {
                    throw std::domain_error(
                        "computeTotalLagrangianElementResponse: material response must be finite");
                }

                // 4. 参考体积积分。materialB 对应 Green-Lagrange 应变的变分算子，
                // 由 referenceBMatrices[a]*F^T 横向拼成；UL 对应的算子是 spatialB。
                Eigen::Matrix<double, 6, Eigen::Dynamic> materialB(6, 3 * nodes);
                for (int a = 0; a < nodes; ++a) {
                    materialB.block<6, 3>(0, 3 * a).noalias() =
                        reference.referenceBMatrices[static_cast<std::size_t>(a)]
                        * kinematics.deformationGradient.transpose();
                }

                // ---------- 材料刚度 K_mat = materialB^T*D*materialB ----------
                // D 对一般超弹性材料会随构形和高斯点改变，也不保证始终正定，
                // 因而在当前积分点直接使用材料返回的相容切线，不再预做常量 LLT。
                const Eigen::MatrixXd constitutiveStiffness =
                    (materialB.transpose() * response.materialTangent * materialB).eval();
                result.tangentStiffness.noalias() +=
                    reference.integrationWeight * constitutiveStiffness;

                // ---------- 几何刚度 K_geo = (G*S*G^T) tensor I3 ----------
                // initialStress 大小为 n_en x n_en；每个标量只加到对应 3x3 块的对角线。
                const Eigen::MatrixXd initialStress =
                    (reference.referenceGradients * response.secondPiolaStress
                     * reference.referenceGradients.transpose()).eval();
                for (int a = 0; a < nodes; ++a) {
                    for (int b = 0; b <= a; ++b) {
                        const double contribution = reference.integrationWeight
                                                    * initialStress(a, b);
                        // &initialStress(a, b)*I加到(a, b)对应的3*3的块上
                        for (int component = 0; component < 3; ++component) {
                            result.tangentStiffness(3 * a + component,
                                                    3 * b + component) += contribution;
                        }
                    }
                }

                // ---------- 内力 fe_A = integral(F*S*G_A)dV ----------
                // F*S 为 3x3，G_A 为 3x1，因此每个局部控制点贡献一个 3x1 向量。
                const Eigen::Matrix3d firstPiola =
                    kinematics.deformationGradient * response.secondPiolaStress;
                for (int a = 0; a < nodes; ++a) {
                    const Eigen::Vector3d gradient =
                        reference.referenceGradients.row(a).transpose();
                    result.internalForce.segment<3>(3 * a).noalias() +=
                        reference.integrationWeight * firstPiola * gradient;
                }
            }
        }
    }
    // 初应力项只累加了下三角；最后显式镜像，得到完整对称 Ke。
    result.tangentStiffness =
        result.tangentStiffness.selfadjointView<Eigen::Lower>();
    return result; // &返回一个element内的切线刚度矩阵+内力向量
}

} // namespace iga
