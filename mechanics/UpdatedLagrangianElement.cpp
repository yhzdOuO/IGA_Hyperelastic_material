#include "mechanics/UpdatedLagrangianElement.h"

#include <Eigen/LU>
#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

ReferenceGradientData evaluateReferenceGradients(
    const nurbs::ShapeFunctionData& shape,
    const nurbs::ControlPointMatrix& localReferencePoints,
    double parameterIntegrationWeight)
{
    const Eigen::Index nodes = localReferencePoints.rows();
    if (nodes == 0 || shape.values.size() != nodes
        || shape.parametricGradients.rows() != nodes
        || !shape.values.allFinite() || !shape.parametricGradients.allFinite()
        || !localReferencePoints.allFinite()
        || !std::isfinite(parameterIntegrationWeight) || parameterIntegrationWeight <= 0.0) {
        throw std::invalid_argument("evaluateReferenceGradients: invalid shape, geometry or weight");
    }
    // ---------- 第一步：参数域 -> 参考构形（书中 9.6--9.8） ----------
    // A0=dX/dxi=X_e^T*dN/dxi。参考映射仅作临时量，UL 不保存 TL 的 B_A。
    const Eigen::Matrix3d referenceJacobian =
        localReferencePoints.transpose() * shape.parametricGradients;
    const double determinant = referenceJacobian.determinant();
    if (!std::isfinite(determinant) || determinant <= 0.0) {
        throw std::domain_error("evaluateReferenceGradients: require det(dX/dxi)>0");
    }
    ReferenceGradientData result;
    result.shapeValues = shape.values;
    result.parametricGradients = shape.parametricGradients;
    // 梯度按行保存，因此右乘 A0^{-1}；按列保存时才是 A0^{-T}。
    result.referenceGradients = shape.parametricGradients * referenceJacobian.inverse();
    result.referenceJacobianDeterminant = determinant;
    result.parameterIntegrationWeight = parameterIntegrationWeight;
    return result;
}

SpatialPointData evaluateSpatialPoint(
    const ReferenceGradientData& reference,
    const nurbs::ControlPointMatrix& localDisplacements,
    const nurbs::ControlPointMatrix& localCurrentPoints)
{
    const Eigen::Index nodes = reference.referenceGradients.rows();
    if (nodes == 0 || reference.parametricGradients.rows() != nodes
        || localDisplacements.rows() != nodes || localCurrentPoints.rows() != nodes
        || !reference.parametricGradients.allFinite() || !localCurrentPoints.allFinite()
        || !std::isfinite(reference.referenceJacobianDeterminant)
        || reference.referenceJacobianDeterminant <= 0.0
        || !std::isfinite(reference.parameterIntegrationWeight)
        || reference.parameterIntegrationWeight <= 0.0) {
        throw std::invalid_argument("evaluateSpatialPoint: inconsistent point data");
    }
    // ---------- 第二步：参数域 -> 当前构形（书中 9.5、9.11） ----------
    SpatialPointData result;
    result.currentJacobian.noalias() = localCurrentPoints.transpose()
                                       * reference.parametricGradients;
    result.currentJacobianDeterminant = result.currentJacobian.determinant();
    if (!std::isfinite(result.currentJacobianDeterminant)
        || result.currentJacobianDeterminant <= 0.0) {
        throw std::domain_error("evaluateSpatialPoint: require det(dx/dxi)>0");
    }
    result.spatialGradients = reference.parametricGradients * result.currentJacobian.inverse();
    // F=I+u_e^T*grad_X(N) 与 x_e^T*grad_X(N) 等价，小位移时更稳健。
    result.kinematics = evaluateKinematics(localDisplacements, reference.referenceGradients);
    result.deformationJacobian = result.kinematics.deformationGradient.determinant();
    if (!std::isfinite(result.deformationJacobian) || result.deformationJacobian <= 0.0) {
        throw std::domain_error("evaluateSpatialPoint: require det(F)>0");
    }
    // 母单元 -> 参数域的因子只乘一次；不能再乘 J，否则会重复计算体积变换。
    result.integrationWeight = result.currentJacobianDeterminant
                               * reference.parameterIntegrationWeight;
    if (!result.spatialGradients.allFinite() || !std::isfinite(result.integrationWeight)) {
        throw std::domain_error("evaluateSpatialPoint: non-finite spatial gradients or weight");
    }
    return result;
}

ElementResult computeUpdatedLagrangianElementResponse(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const nurbs::ControlPointMatrix& localDisplacements,
    const HyperelasticMaterial& material)
{
    const int nodes = space.numElementNodes();
    if (localDisplacements.rows() != nodes || !localDisplacements.allFinite()) {
        throw std::invalid_argument(
            "computeUpdatedLagrangianElementResponse: invalid local displacements");
    }
    // 控制点收集和 x_A=X_A+u_A 在积分循环外完成。
    const nurbs::ControlPointMatrix referencePoints = space.elementControlPoints(element);
    const nurbs::ControlPointMatrix currentPoints = referencePoints + localDisplacements;
    const double parentJacobian = space.parentJacobianDeterminant(element);

    ElementResult result;
    result.tangentStiffness = Eigen::MatrixXd::Zero(3 * nodes, 3 * nodes);
    // 仅初始化积分累加器；实际内力在下方 Gauss 循环中按式 (9.15b) 累加。
    result.internalForce = Eigen::VectorXd::Zero(3 * nodes);
    Eigen::Matrix<double, 6, Eigen::Dynamic> spatialB(6, 3 * nodes);
    for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
        for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
            for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                const std::array<int, 3> gaussIndices{q0, q1, q2};
                // 1. N、参数梯度和参考梯度。
                const ReferenceGradientData reference = evaluateReferenceGradients(
                    space.evaluateShapeFunctions(element, gaussIndices), referencePoints,
                    parentJacobian * space.gaussWeight(gaussIndices));
                // 2. 当前梯度、当前积分权重，以及总 F、J。
                const SpatialPointData point = evaluateSpatialPoint(
                    reference, localDisplacements, currentPoints);
                // 3. 材料提供 Cauchy 应力和空间弹性矩阵。
                const SpatialMaterialResponse response = material.evaluateSpatial(point.kinematics);
                if (!response.cauchyStress.allFinite() || !response.spatialTangent.allFinite()
                    || !std::isfinite(response.strainEnergyDensity)) {
                    throw std::domain_error(
                        "computeUpdatedLagrangianElementResponse: non-finite material response");
                }

                // 4. 当前体积积分（书中 9.15b、9.41、9.44）。
                for (int a = 0; a < nodes; ++a) {
                    const Eigen::Vector3d gradient = point.spatialGradients.row(a).transpose();
                    spatialB.block<6, 3>(0, 3 * a) = makeLinearB(gradient);
                    // 式 (9.15b)：f_int,a = integral_{v_e}(sigma * grad_x(N_a)) dv。
                    // 每个积分点累加 w_q*sigma_q*grad_x(N_a)，全局层再形成 R=f_int-f_ext。
                    result.internalForce.segment<3>(3 * a).noalias() +=
                        point.integrationWeight * response.cauchyStress * gradient;
                }
                // 切线刚度
                result.tangentStiffness.noalias() += point.integrationWeight
                    * spatialB.transpose() * response.spatialTangent * spatialB;
                const Eigen::MatrixXd initialStress = point.spatialGradients
                    * response.cauchyStress * point.spatialGradients.transpose();
                for (int a = 0; a < nodes; ++a) {
                    for (int b = 0; b < nodes; ++b) {
                        result.tangentStiffness.block<3, 3>(3 * a, 3 * b).diagonal().array() +=
                            point.integrationWeight * initialStress(a, b);
                    }
                }
            }
        }
    }
    return result;
}

} // namespace iga
