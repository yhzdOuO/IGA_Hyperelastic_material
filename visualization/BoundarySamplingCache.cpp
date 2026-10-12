#include "visualization/BoundarySamplingCache.h"

#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <Eigen/Eigenvalues>
#include <cmath>
#include <stdexcept>

namespace visualization {

BoundarySamplingCache::BoundarySamplingCache(
    const nurbs::TensorProductNurbsSpace& space, int samplesPerSpan)
    : dofCount_(3 * static_cast<int>(space.controlPoints().rows()))
{
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(dofCount_);
    const iga::SaintVenantKirchhoffMaterial dummyMaterial(0.0, 1.0);
    const FieldEvaluator evaluator(space, zero, dummyMaterial,
                                   Eigen::Vector3d::UnitX(),
                                   -Eigen::Vector3d::UnitZ());
    referenceMesh_ = BoundarySampler(samplesPerSpan).sample(space, evaluator);
    operators_.reserve(referenceMesh_.vertices.size());
    for (const SurfaceVertex& vertex : referenceMesh_.vertices) {
        const nurbs::PointShapeFunctionData basis =
            space.evaluateShapeFunctionsAt(vertex.parameters);
        SampleOperator sample;
        sample.connectivity = basis.connectivity;
        sample.values = basis.shape.values;
        sample.parametricGradients = basis.shape.parametricGradients;
        nurbs::ControlPointMatrix points(sample.connectivity.size(), 3);
        for (std::size_t local = 0; local < sample.connectivity.size(); ++local) {
            points.row(static_cast<int>(local)) =
                space.controlPoints().row(sample.connectivity[local]);
        }
        const Eigen::Matrix3d jacobian = points.transpose()
            * sample.parametricGradients;
        if (!std::isfinite(jacobian.determinant()) || jacobian.determinant() <= 0.0) {
            throw std::domain_error("BoundarySamplingCache: invalid reference Jacobian");
        }
        sample.physicalGradients = sample.parametricGradients * jacobian.inverse();
        operators_.push_back(std::move(sample));
    }
}

std::vector<SampleScalars> BoundarySamplingCache::evaluateScalars(
    const Eigen::VectorXd& displacements,
    const iga::HyperelasticMaterial& material,
    const Eigen::Vector3d& axialDirection,
    const Eigen::Vector3d& gravityDirection) const
{
    if (displacements.size() != dofCount_ || !displacements.allFinite()) {
        throw std::invalid_argument("BoundarySamplingCache: invalid displacement vector");
    }
    const Eigen::Vector3d axial = axialDirection.normalized();
    const Eigen::Vector3d gravity = gravityDirection.normalized();
    std::vector<SampleScalars> result(operators_.size());
    for (std::size_t index = 0; index < operators_.size(); ++index) {
        const SampleOperator& sample = operators_[index];
        nurbs::ControlPointMatrix local(sample.connectivity.size(), 3);
        for (std::size_t i = 0; i < sample.connectivity.size(); ++i) {
            local.row(static_cast<int>(i)) =
                displacements.segment<3>(3 * sample.connectivity[i]).transpose();
        }
        // 功能块 1：缓存梯度只负责生成 H、F、E，材料模型随后给出第二 PK 应力。
        const iga::KinematicState state = iga::evaluateKinematics(
            local, sample.physicalGradients);
        const iga::MaterialResponse response = material.evaluate(state);
        const double determinant = state.deformationGradient.determinant();
        if (!std::isfinite(determinant) || determinant <= 0.0) {
            throw std::domain_error("BoundarySamplingCache: invalid deformation Jacobian");
        }
        const Eigen::Matrix3d sigma = state.deformationGradient
            * response.secondPiolaStress * state.deformationGradient.transpose()
            / determinant;
        SampleScalars& values = result[index];
        values[static_cast<int>(ScalarField::AxialStress)] = axial.dot(sigma * axial);
        values[static_cast<int>(ScalarField::GravityStress)] = gravity.dot(sigma * gravity);
        values[static_cast<int>(ScalarField::ShearXY)] = sigma(0, 1);
        values[static_cast<int>(ScalarField::ShearYZ)] = sigma(1, 2);
        values[static_cast<int>(ScalarField::ShearXZ)] = sigma(0, 2);
        const Eigen::Matrix3d deviator = sigma
            - (sigma.trace() / 3.0) * Eigen::Matrix3d::Identity();
        values[static_cast<int>(ScalarField::VonMises)] =
            std::sqrt(1.5 * deviator.squaredNorm());
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(sigma);
        if (eigenSolver.info() != Eigen::Success) {
            throw std::runtime_error("BoundarySamplingCache: stress eigensolve failed");
        }
        values[static_cast<int>(ScalarField::PrincipalMin)] = eigenSolver.eigenvalues()(0);
        values[static_cast<int>(ScalarField::PrincipalMax)] = eigenSolver.eigenvalues()(2);
    }
    return result;
}

SurfaceMesh BoundarySamplingCache::makeFrame(
    const Eigen::VectorXd& displacements,
    const std::vector<SampleScalars>& scalars) const
{
    if (displacements.size() != dofCount_ || !displacements.allFinite()
        || scalars.size() != operators_.size()) {
        throw std::invalid_argument("BoundarySamplingCache: invalid frame data");
    }
    SurfaceMesh mesh = referenceMesh_;
    for (std::size_t index = 0; index < operators_.size(); ++index) {
        const SampleOperator& sample = operators_[index];
        nurbs::ControlPointMatrix local(sample.connectivity.size(), 3);
        for (std::size_t i = 0; i < sample.connectivity.size(); ++i) {
            local.row(static_cast<int>(i)) =
                displacements.segment<3>(3 * sample.connectivity[i]).transpose();
        }
        SurfaceVertex& vertex = mesh.vertices[index];
        vertex.displacement = local.transpose() * sample.values;
        const Eigen::Matrix3d gradient = local.transpose()
            * sample.parametricGradients;
        for (int direction = 0; direction < 3; ++direction) {
            vertex.displacementTangents[direction] = gradient.col(direction);
        }
        vertex.scalars = scalars[index];
    }
    mesh.updateGeometry(1.0);
    mesh.selectField(ScalarField::AxialStress);
    return mesh;
}

} // namespace visualization
