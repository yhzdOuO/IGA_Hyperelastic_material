#include "mechanics/TotalLagrangianElement.h"

#include <cmath>
#include <stdexcept>

namespace iga {

BMatrix makeLinearB(const Eigen::Vector3d& physicalGradient)
{
    // G_A 是 3x1 的物理域梯度；B_A 是 6x3，前三行对应正应变，
    // 后三行对应工程剪应变 2E12、2E23、2E13。
    const double g1 = physicalGradient(0);
    const double g2 = physicalGradient(1);
    const double g3 = physicalGradient(2);
    BMatrix result;
    result << g1, 0.0, 0.0,
              0.0, g2, 0.0,
              0.0, 0.0, g3,
              g2, g1, 0.0,
              0.0, g3, g2,
              g3, 0.0, g1;
    return result;
}

ReferencePointData evaluateReferencePoint(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const std::array<int, 3>& gaussIndices)
{
    // ---------- L1：参考构形中的几何量（Newton 迭代中保持不变） ----------
    const nurbs::ShapeFunctionData shape =
        space.evaluateShapeFunctions(element, gaussIndices); // &计算型函数和参数梯度: R 和 \nabla R
    const nurbs::ControlPointMatrix localControlPoints =
        space.elementControlPoints(element);

    ReferencePointData result;
    result.shapeValues = shape.values;
    result.parametricGradients = shape.parametricGradients;
    // P_e: n_en x 3，dN/dxi: n_en x 3，因此
    // J = P_e^T * dN/dxi 的大小为 3x3，J_ij=dX_i/dxi_j。
    result.jacobian.noalias() = localControlPoints.transpose()
                                * result.parametricGradients;
    result.jacobianDeterminant = result.jacobian.determinant();
    if (!std::isfinite(result.jacobianDeterminant)
        || result.jacobianDeterminant <= 0.0) {
        throw std::domain_error(
            "evaluateReferencePoint: geometry Jacobian determinant must be positive");
    }
    result.inverseJacobian = result.jacobian.inverse();
    // dN/dX = dN/dxi * J^{-1}，结果仍为 n_en x 3；注意这里不转置。
    result.physicalGradients.noalias() = result.parametricGradients
                                         * result.inverseJacobian;
    // 积分权重同时包含几何 Jacobian、母单元 Jacobian 和三维 Gauss 权重。
    result.integrationWeight = result.jacobianDeterminant
        * space.parentJacobianDeterminant(element)
        * space.gaussWeight(gaussIndices); //& total lagrangian下，当前构型改变不用更新Jacobian，该Jacobian是参数空间到参考构型的映射，并非到当前构型

    // 每个局部基函数保存一个 6x3 的 B_A，后续各 Newton 迭代可复用。
    result.linearBMatrices.reserve(static_cast<std::size_t>(space.numElementNodes()));
    for (int a = 0; a < space.numElementNodes(); ++a) {
        const Eigen::Vector3d gradient = result.physicalGradients.row(a).transpose();
        result.linearBMatrices.push_back(makeLinearB(gradient));
    }
    return result; //型函数+参数梯度，参数空间到物理空间的Jacobian矩阵及其行列式，物理梯度，积分权重
}

KinematicState evaluateKinematics(
    const nurbs::ControlPointMatrix& localDisplacements,
    const GradientMatrix& physicalGradients)
{
    if (localDisplacements.rows() != physicalGradients.rows()
        || localDisplacements.cols() != 3 || physicalGradients.cols() != 3
        || !localDisplacements.allFinite() || !physicalGradients.allFinite()) {
        throw std::invalid_argument("evaluateKinematics: displacement dimensions are inconsistent");
    }
    // ---------- L2：随当前位移变化的运动学量 ----------
    KinematicState result;
    // u_e: n_en x 3，G=dN/dX: n_en x 3，故 H=u_e^T*G 为 3x3。
    result.displacementGradient.noalias() = localDisplacements.transpose()
                                            * physicalGradients;
    result.deformationGradient = Eigen::Matrix3d::Identity()
                                 + result.displacementGradient;
    // 用 H 形式计算 E，避免小位移下 F^T*F-I 的浮点抵消。
    const Eigen::Matrix3d quadratic = result.displacementGradient.transpose()
                                      * result.displacementGradient;
    result.greenLagrangeStrain = 0.5 * (result.displacementGradient
                                        + result.displacementGradient.transpose()
                                        + quadratic);
    return result;
}

Eigen::Matrix<double, 6, 1> stressVoigt(const Eigen::Matrix3d& stress)
{
    // 应力 Voigt 向量为 6x1，剪切分量不乘 2，以保证 S^T*deltaE = S:deltaE。
    Eigen::Matrix<double, 6, 1> result;
    result << stress(0, 0), stress(1, 1), stress(2, 2),
              stress(0, 1), stress(1, 2), stress(0, 2);
    return result;
}

ElementResult computeElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material)
{
    const int nodes = space.numElementNodes();
    if (localDisplacements.rows() != nodes || localDisplacements.cols() != 3
        || !localDisplacements.allFinite()) {
        throw std::invalid_argument(
            "computeElementResponse: local displacement dimensions are inconsistent");
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
                const ReferencePointData reference =
                    evaluateReferencePoint(space, element, gaussIndices); // 计算型函数，NURBS相关
                const KinematicState state = evaluateKinematics(
                    localDisplacements, reference.physicalGradients); // 位移d相关的运动学量
                // 运动学和材料计算分成两步；这样同一套单元积分可以复用任意超弹性模型。
                const MaterialResponse response = material.evaluate(state);
                if (!response.secondPiolaStress.allFinite()
                    || !response.materialTangent.allFinite()
                    || !std::isfinite(response.strainEnergyDensity)) {
                    throw std::domain_error(
                        "computeElementResponse: material response must be finite");
                }

                // 将每个 Btilde_A=B_A*F^T（6x3）横向拼成
                // Btilde（6 x 3*n_en），每个控制点只计算一次。
                Eigen::Matrix<double, 6, Eigen::Dynamic> nonlinearB(6, 3 * nodes);
                for (int a = 0; a < nodes; ++a) {
                    nonlinearB.block<6, 3>(0, 3 * a).noalias() =
                        reference.linearBMatrices[static_cast<std::size_t>(a)]
                        * state.deformationGradient.transpose();
                }

                // ---------- 材料刚度 K_mat = Btilde^T*D*Btilde ----------
                // D 对一般超弹性材料会随构形和高斯点改变，也不保证始终正定，
                // 因而在当前积分点直接使用材料返回的相容切线，不再预做常量 LLT。
                const Eigen::MatrixXd materialTangent =
                    (nonlinearB.transpose() * response.materialTangent * nonlinearB).eval();
                result.tangentStiffness.noalias() +=
                    reference.integrationWeight * materialTangent;

                // ---------- 几何刚度 K_geo = (G*S*G^T) tensor I3 ----------
                // geometricScalar 大小为 n_en x n_en；每个标量只加到对应 3x3 块的对角线。
                const Eigen::MatrixXd geometricScalar =
                    (reference.physicalGradients * response.secondPiolaStress
                     * reference.physicalGradients.transpose()).eval();
                for (int a = 0; a < nodes; ++a) {
                    for (int b = 0; b <= a; ++b) {
                        const double contribution = reference.integrationWeight
                                                    * geometricScalar(a, b);
                        // &geometricScalar(a, b)*I加到(a, b)对应的3*3的块上
                        for (int component = 0; component < 3; ++component) {
                            result.tangentStiffness(3 * a + component,
                                                    3 * b + component) += contribution;
                        }
                    }
                }

                // ---------- 内力 fe_A = integral(F*S*G_A)dV ----------
                // F*S 为 3x3，G_A 为 3x1，因此每个局部控制点贡献一个 3x1 向量。
                const Eigen::Matrix3d firstPiola =
                    state.deformationGradient * response.secondPiolaStress;
                for (int a = 0; a < nodes; ++a) {
                    const Eigen::Vector3d gradient =
                        reference.physicalGradients.row(a).transpose();
                    result.internalForce.segment<3>(3 * a).noalias() +=
                        reference.integrationWeight * firstPiola * gradient;
                }
            }
        }
    }
    // rankUpdate 和几何项只写了下三角；最后显式镜像，得到完整对称 Ke。
    result.tangentStiffness =
        result.tangentStiffness.selfadjointView<Eigen::Lower>();
    return result; // &返回一个element内的切线刚度矩阵+内力向量
}

} // namespace iga
