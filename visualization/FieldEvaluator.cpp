#include "visualization/FieldEvaluator.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>

namespace visualization {

FieldEvaluator::FieldEvaluator(const nurbs::TensorProductNurbsSpace& space,
                               const Eigen::VectorXd& globalDisplacements,
                               const iga::HyperelasticMaterial& material,
                               const Eigen::Vector3d& axialDirection,
                               const Eigen::Vector3d& gravityDirection)
    : space_(space),
      displacements_(globalDisplacements),
      material_(material),
      axialDirection_(axialDirection),
      gravityDirection_(gravityDirection)
{
    if (displacements_.size() != 3 * space_.controlPoints().rows()
        || !displacements_.allFinite()) {
        throw std::invalid_argument("FieldEvaluator: displacement vector has the wrong size");
    }
    if (!axialDirection_.allFinite() || axialDirection_.norm() <= 0.0
        || !gravityDirection_.allFinite() || gravityDirection_.norm() <= 0.0) {
        throw std::invalid_argument("FieldEvaluator: directions must be finite and nonzero");
    }
    axialDirection_.normalize();
    gravityDirection_.normalize();
}

FieldPoint FieldEvaluator::evaluate(const std::array<double, 3>& parameters) const
{
    // 功能块 1：复用 NURBS 空间的任意点基函数与参数导数，禁止另写一套算法。
    const nurbs::PointShapeFunctionData basis =
        space_.evaluateShapeFunctionsAt(parameters);
    const int count = static_cast<int>(basis.connectivity.size());
    nurbs::ControlPointMatrix points(count, 3);
    nurbs::ControlPointMatrix displacement(count, 3);
    // 根据活动基函数连接关系，收集局部控制点坐标和真实控制点位移。
    for (int local = 0; local < count; ++local) {
        const int global = basis.connectivity[static_cast<std::size_t>(local)];
        points.row(local) = space_.controlPoints().row(global);
        displacement.row(local) = displacements_.segment<3>(3 * global).transpose();
    }

    // 功能块 2：插值参考位置、位移，并形成 dX/dxi 与 du/dxi。
    FieldPoint result;
    result.referencePosition = points.transpose() * basis.shape.values;
    result.displacement = displacement.transpose() * basis.shape.values;
    const Eigen::Matrix3d referenceJacobian =
        points.transpose() * basis.shape.parametricGradients;
    // 正 Jacobian 是后续求物理梯度及应力的必要条件。
    const double referenceDeterminant = referenceJacobian.determinant();
    if (!std::isfinite(referenceDeterminant) || referenceDeterminant <= 0.0) {
        throw std::domain_error("FieldEvaluator: reference Jacobian must be positive");
    }
    const Eigen::Matrix3d displacementParametricGradient =
        displacement.transpose() * basis.shape.parametricGradients;
    // 保存三个参数方向的切向，供不同边界面的解析法向和形变放大使用。
    for (int direction = 0; direction < 3; ++direction) {
        result.referenceTangents[direction] = referenceJacobian.col(direction);
        result.displacementTangents[direction] =
            displacementParametricGradient.col(direction);
    }

    // 功能块 3：商法则后的参数导数通过 J^{-1} 转成参考物理坐标梯度。
    const iga::GradientMatrix physicalGradients =
        basis.shape.parametricGradients * referenceJacobian.inverse();
    // 功能块 3：先复用求解器得到 H、F、E，再由通用材料接口单独计算 S。
    const iga::KinematicState state =
        iga::evaluateKinematics(displacement, physicalGradients);
    const iga::MaterialResponse response = material_.evaluate(state);
    const double deformationJacobian = state.deformationGradient.determinant();
    if (!std::isfinite(deformationJacobian) || deformationJacobian <= 0.0) {
        throw std::domain_error("FieldEvaluator: deformation Jacobian must be positive");
    }
    // 关键步骤：Total Lagrangian 求得的是 PK2，应推前为 sigma=F*S*F^T/det(F)。
    result.cauchyStress = (state.deformationGradient * response.secondPiolaStress
                           * state.deformationGradient.transpose())
                          / deformationJacobian;

    // 功能块 4：一次计算并缓存所有可切换标量，切换字段时无需重新采样。
    const Eigen::Matrix3d& sigma = result.cauchyStress;
    result.scalars[static_cast<int>(ScalarField::AxialStress)] =
        axialDirection_.dot(sigma * axialDirection_);
    result.scalars[static_cast<int>(ScalarField::GravityStress)] =
        gravityDirection_.dot(sigma * gravityDirection_);
    result.scalars[static_cast<int>(ScalarField::ShearXY)] = sigma(0, 1);
    result.scalars[static_cast<int>(ScalarField::ShearYZ)] = sigma(1, 2);
    result.scalars[static_cast<int>(ScalarField::ShearXZ)] = sigma(0, 2);
    // von Mises 使用 Cauchy 应力偏量的三维定义。
    const double mean = sigma.trace() / 3.0;
    const Eigen::Matrix3d deviator = sigma - mean * Eigen::Matrix3d::Identity();
    result.scalars[static_cast<int>(ScalarField::VonMises)] =
        std::sqrt(1.5 * deviator.squaredNorm());
    // Cauchy 应力为实对称矩阵，自伴特征求解器给出升序主应力。
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(sigma);
    if (eigenSolver.info() != Eigen::Success) {
        throw std::runtime_error("FieldEvaluator: principal stress calculation failed");
    }
    result.scalars[static_cast<int>(ScalarField::PrincipalMin)] =
        eigenSolver.eigenvalues()(0);
    result.scalars[static_cast<int>(ScalarField::PrincipalMax)] =
        eigenSolver.eigenvalues()(2);
    return result;
}

} // namespace visualization
