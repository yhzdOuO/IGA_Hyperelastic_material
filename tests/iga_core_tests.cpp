#include "mechanics/TotalLagrangianElement.h"
#include "mechanics/SaintVenantKirchhoffMaterial.h"
#include "nurbs/NurbsAlgorithm.h"

#include <Eigen/Geometry>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using nurbs::ControlPointMatrix;
using nurbs::KnotVector;
using nurbs::TensorProductNurbsSpace;

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void requireNear(double actual, double expected, double tolerance,
                 const std::string& message)
{
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": got " + std::to_string(actual)
                                 + ", expected " + std::to_string(expected));
    }
}

TensorProductNurbsSpace makeUnitCube(int degree, bool rational = false)
{
    std::vector<double> knots(static_cast<std::size_t>(2 * degree + 2));
    for (int i = 0; i <= degree; ++i) {
        knots[static_cast<std::size_t>(i)] = 0.0;
        knots[static_cast<std::size_t>(degree + 1 + i)] = 1.0;
    }
    const int count = degree + 1;
    const int total = count * count * count;
    ControlPointMatrix points(total, 3);
    Eigen::VectorXd weights = Eigen::VectorXd::Ones(total);
    int global = 0;
    for (int k = 0; k < count; ++k) {
        for (int j = 0; j < count; ++j) {
            for (int i = 0; i < count; ++i) {
                points.row(global) << static_cast<double>(i) / degree,
                                      static_cast<double>(j) / degree,
                                      static_cast<double>(k) / degree;
                if (rational) {
                    weights(global) = 0.85 + 0.03 * global;
                }
                ++global;
            }
        }
    }
    return TensorProductNurbsSpace(
        {KnotVector(knots), KnotVector(knots), KnotVector(knots)},
        {degree, degree, degree}, std::move(points), std::move(weights));
}

TensorProductNurbsSpace makeTwoByTwoByTwoLinearSpace()
{
    const std::vector<double> knots{0.0, 0.0, 0.5, 1.0, 1.0};
    ControlPointMatrix points(27, 3);
    int global = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                points.row(global++) << 0.5 * i, 0.5 * j, 0.5 * k;
            }
        }
    }
    return TensorProductNurbsSpace(
        {KnotVector(knots), KnotVector(knots), KnotVector(knots)},
        {1, 1, 1}, std::move(points), Eigen::VectorXd::Ones(27));
}

void testDerivativeStorageAndValues()
{
    const KnotVector knots({0.0, 0.0, 0.0, 1.0, 1.0, 1.0});
    Eigen::MatrixXd derivatives;
    nurbs::DersBasisFuns(2, 0.5, 2, 3, knots, derivatives);
    require(derivatives.rows() == 4 && derivatives.cols() == 3,
            "DersBasisFuns must return a contiguous Eigen matrix");
    requireNear(derivatives.row(0).sum(), 1.0, 1e-14, "basis partition of unity");
    requireNear(derivatives.row(1).sum(), 0.0, 1e-14, "first derivatives sum to zero");
    require(derivatives.row(3).isZero(0.0), "derivatives above degree must be zero");
}

void testSpansAndLazyIEN()
{
    const TensorProductNurbsSpace space = makeTwoByTwoByTwoLinearSpace();
    require(space.numElements() == 8, "nonzero spans must form 8 tensor-product elements");
    require(space.span(0, 0) == 1 && space.span(0, 1) == 2,
            "directional span table is incorrect");
    const std::vector<int> expectedFirst{0, 1, 3, 4, 9, 10, 12, 13};
    const std::vector<int> expectedLast{13, 14, 16, 17, 22, 23, 25, 26};
    require(space.elementConnectivity(0) == expectedFirst, "first element IEN is incorrect");
    require(space.elementConnectivity(7) == expectedLast, "last element IEN is incorrect");
    for (int element = 0; element < space.numElements(); ++element) {
        const std::vector<int> connectivity = space.elementConnectivity(element);
        for (int local = 0; local < space.numElementNodes(); ++local) {
            require(space.IEN(local, element) == connectivity[static_cast<std::size_t>(local)],
                    "lazy IEN(a,e) differs from connectivity view");
        }
    }
}

void testShapeFunctionsAndRationalDerivatives()
{
    const TensorProductNurbsSpace polynomial = makeUnitCube(2);
    const nurbs::ShapeFunctionData center =
        polynomial.evaluateShapeFunctions(0, {1, 1, 1});
    const double oneDimensional[3] = {0.25, 0.5, 0.25};
    int local = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                requireNear(center.values(local++),
                            oneDimensional[i] * oneDimensional[j] * oneDimensional[k],
                            2e-14, "unit weights must recover tensor-product B-splines");
            }
        }
    }
    requireNear(center.values.sum(), 1.0, 2e-14, "NURBS values must sum to one");
    require(center.parametricGradients.colwise().sum().norm() < 2e-14,
            "NURBS parametric gradients must sum to zero");

    const TensorProductNurbsSpace rational = makeUnitCube(2, true);
    const nurbs::ShapeFunctionData rationalCenter =
        rational.evaluateShapeFunctions(0, {1, 1, 1});
    requireNear(rationalCenter.values.sum(), 1.0, 2e-14,
                "rational values must sum to one");
    require(rationalCenter.parametricGradients.colwise().sum().norm() < 3e-14,
            "rational quotient-rule gradients must sum to zero");

    const double derivative[3] = {-1.0, 0.0, 1.0};
    Eigen::VectorXd weighted(27);
    Eigen::Matrix<double, Eigen::Dynamic, 3> weightedDerivative(27, 3);
    double sum = 0.0;
    Eigen::Vector3d derivativeSum = Eigen::Vector3d::Zero();
    local = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                const double weight = rational.weights()(local);
                weighted(local) = weight * oneDimensional[i] * oneDimensional[j]
                                  * oneDimensional[k];
                weightedDerivative.row(local) <<
                    weight * derivative[i] * oneDimensional[j] * oneDimensional[k],
                    weight * oneDimensional[i] * derivative[j] * oneDimensional[k],
                    weight * oneDimensional[i] * oneDimensional[j] * derivative[k];
                sum += weighted(local);
                derivativeSum += weightedDerivative.row(local).transpose();
                ++local;
            }
        }
    }
    for (int a = 0; a < 27; ++a) {
        const double expectedValue = weighted(a) / sum;
        const Eigen::Vector3d expectedGradient =
            (weightedDerivative.row(a).transpose() - expectedValue * derivativeSum) / sum;
        requireNear(rationalCenter.values(a), expectedValue, 3e-14,
                    "rational quotient-rule value is incorrect");
        require((rationalCenter.parametricGradients.row(a).transpose()
                 - expectedGradient).norm() < 5e-14,
                "rational quotient-rule gradient is incorrect");
    }
}

void testJacobiansAndPhysicalGradients()
{
    const TensorProductNurbsSpace space = makeUnitCube(2, true);
    for (int q2 = 0; q2 < 3; ++q2) {
        for (int q1 = 0; q1 < 3; ++q1) {
            for (int q0 = 0; q0 < 3; ++q0) {
                const iga::ReferencePointData point =
                    iga::evaluateReferencePoint(space, 0, {q0, q1, q2});
                require(point.jacobianDeterminant > 0.0, "det(J) must be positive");
                const Eigen::Matrix3d identity =
                    space.elementControlPoints(0).transpose() * point.physicalGradients;
                require((identity - Eigen::Matrix3d::Identity()).norm() < 3e-13,
                        "P^T grad_X(N) must equal I");
                require(point.physicalGradients.colwise().sum().norm() < 2e-13,
                        "physical gradients must sum to zero");
            }
        }
    }
    const TensorProductNurbsSpace refined = makeTwoByTwoByTwoLinearSpace();
    const iga::ReferencePointData point = iga::evaluateReferencePoint(refined, 0, {0, 0, 0});
    requireNear(refined.parentJacobianDeterminant(0), 1.0 / 64.0, 1e-15,
                "parent-element Jacobian must be included explicitly");
    requireNear(point.integrationWeight, 1.0 / 64.0, 2e-15,
                "integration weight must contain both Jacobians");
}

void testKinematicsAndBConventions()
{
    const TensorProductNurbsSpace space = makeUnitCube(2);
    const iga::ReferencePointData point = iga::evaluateReferencePoint(space, 0, {1, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(7.0, 3.0);
    const Eigen::Matrix3d rotation =
        Eigen::AngleAxisd(0.37, Eigen::Vector3d(1.0, 2.0, -1.0).normalized())
            .toRotationMatrix();
    const ControlPointMatrix rigidDisplacement =
        (space.elementControlPoints(0)
         * (rotation.transpose() - Eigen::Matrix3d::Identity())).eval();
    const iga::KinematicState rigid = iga::evaluateKinematics(
        rigidDisplacement, point.physicalGradients);
    const iga::MaterialResponse rigidResponse = material.evaluate(rigid);
    require((rigid.deformationGradient - rotation).norm() < 3e-14,
            "rigid rotation must produce the prescribed F");
    require(rigid.greenLagrangeStrain.norm() < 4e-14,
            "rigid rotation must produce zero Green-Lagrange strain");
    require(rigidResponse.secondPiolaStress.norm() < 6e-13,
            "rigid rotation must produce zero material stress");
    require(std::abs(rigidResponse.strainEnergyDensity) < 1e-26,
            "rigid rotation must produce zero material energy");

    const Eigen::Vector3d gradient(0.4, -0.2, 0.7);
    const iga::BMatrix linearB = iga::makeLinearB(gradient);
    require((linearB * Eigen::Matrix3d::Identity() - linearB).norm() == 0.0,
            "Btilde must reduce to B when F=I");

    ControlPointMatrix displacement = 0.04 * space.elementControlPoints(0);
    displacement.col(1) += 0.03 * space.elementControlPoints(0).col(0);
    const iga::KinematicState state = iga::evaluateKinematics(
        displacement, point.physicalGradients);
    const iga::MaterialResponse materialResponse = material.evaluate(state);
    const Eigen::Matrix<double, 6, 1> stressVector =
        iga::stressVoigt(materialResponse.secondPiolaStress);
    for (int a = 0; a < space.numElementNodes(); ++a) {
        const Eigen::Vector3d g = point.physicalGradients.row(a).transpose();
        const iga::BMatrix nonlinearB =
            point.linearBMatrices[static_cast<std::size_t>(a)]
            * state.deformationGradient.transpose();
        const Eigen::Vector3d voigtForce = nonlinearB.transpose() * stressVector;
        const Eigen::Vector3d tensorForce =
            state.deformationGradient * materialResponse.secondPiolaStress * g;
        require((voigtForce - tensorForce).norm() < 2e-14,
                "Voigt and tensor internal-force paths must agree");
    }
}

Eigen::VectorXd flattenDisplacements(const ControlPointMatrix& displacement)
{
    Eigen::VectorXd result(3 * displacement.rows());
    for (Eigen::Index a = 0; a < displacement.rows(); ++a) {
        result.segment<3>(3 * a) = displacement.row(a).transpose();
    }
    return result;
}

ControlPointMatrix unflattenDisplacements(const Eigen::VectorXd& displacement)
{
    ControlPointMatrix result(displacement.size() / 3, 3);
    for (Eigen::Index a = 0; a < result.rows(); ++a) {
        result.row(a) = displacement.segment<3>(3 * a).transpose();
    }
    return result;
}

void testElementTangentByFiniteDifference()
{
    const TensorProductNurbsSpace space = makeUnitCube(1);
    const iga::SaintVenantKirchhoffMaterial material(5.0, 2.0);
    ControlPointMatrix displacement = ControlPointMatrix::Zero(8, 3);
    const ControlPointMatrix points = space.elementControlPoints(0);
    displacement.col(0) = 0.03 * points.col(0) + 0.01 * points.col(1);
    displacement.col(1) = -0.02 * points.col(1) + 0.015 * points.col(2);
    displacement.col(2) = 0.01 * points.col(0) + 0.025 * points.col(2);
    const iga::ElementResult response =
        iga::computeElementResponse(space, 0, displacement, material);
    require((response.tangentStiffness - response.tangentStiffness.transpose()).norm()
                < 2e-13,
            "element tangent must be symmetric");

    Eigen::VectorXd direction(24);
    for (int i = 0; i < direction.size(); ++i) {
        direction(i) = std::sin(0.7 * (i + 1));
    }
    direction.normalize();
    const Eigen::VectorXd base = flattenDisplacements(displacement);
    const Eigen::VectorXd exact = response.tangentStiffness * direction;
    double previousError = 0.0;
    for (int level = 0; level < 3; ++level) {
        const double step = 2.0e-3 / std::pow(2.0, level);
        const ControlPointMatrix plus = unflattenDisplacements(base + step * direction);
        const ControlPointMatrix minus = unflattenDisplacements(base - step * direction);
        const Eigen::VectorXd difference =
            (iga::computeElementResponse(space, 0, plus, material).internalForce
             - iga::computeElementResponse(space, 0, minus, material).internalForce)
            / (2.0 * step);
        const double error = (difference - exact).norm();
        if (level > 0) {
            require(previousError / error > 3.8,
                    "central-difference tangent check must converge at second order");
        }
        previousError = error;
    }
    require(previousError < 2e-7, "finite-difference tangent error is too large");
}

void testElementResponseAgainstDirectDenseTangent()
{
    const TensorProductNurbsSpace space = makeUnitCube(1);
    const iga::SaintVenantKirchhoffMaterial material(4.0, 1.5);
    ControlPointMatrix displacement = 0.02 * space.elementControlPoints(0);
    const iga::ElementResult optimized =
        iga::computeElementResponse(space, 0, displacement, material);
    Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(24, 24);
    for (int q2 = 0; q2 < 2; ++q2) {
        for (int q1 = 0; q1 < 2; ++q1) {
            for (int q0 = 0; q0 < 2; ++q0) {
                const iga::ReferencePointData reference =
                    iga::evaluateReferencePoint(space, 0, {q0, q1, q2});
                const iga::KinematicState state = iga::evaluateKinematics(
                    displacement, reference.physicalGradients);
                const iga::MaterialResponse materialResponse = material.evaluate(state);
                Eigen::Matrix<double, 6, Eigen::Dynamic> nonlinearB(6, 24);
                for (int a = 0; a < 8; ++a) {
                    nonlinearB.block<6, 3>(0, 3 * a) =
                        reference.linearBMatrices[static_cast<std::size_t>(a)]
                        * state.deformationGradient.transpose();
                }
                dense.noalias() += reference.integrationWeight
                    * nonlinearB.transpose() * materialResponse.materialTangent * nonlinearB;
                const Eigen::MatrixXd geometric =
                    (reference.physicalGradients * materialResponse.secondPiolaStress
                     * reference.physicalGradients.transpose()).eval();
                for (int a = 0; a < 8; ++a) {
                    for (int b = 0; b < 8; ++b) {
                        for (int component = 0; component < 3; ++component) {
                            dense(3 * a + component, 3 * b + component) +=
                                reference.integrationWeight * geometric(a, b);
                        }
                    }
                }
            }
        }
    }
    require((optimized.tangentStiffness - dense).norm() < 2e-13,
            "element and direct dense tangent routes must agree");
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"Eigen derivative storage", testDerivativeStorageAndValues},
        {"span tables and lazy IEN", testSpansAndLazyIEN},
        {"NURBS values and derivatives", testShapeFunctionsAndRationalDerivatives},
        {"two Jacobians and physical gradients", testJacobiansAndPhysicalGradients},
        {"kinematics and Voigt conventions", testKinematicsAndBConventions},
        {"finite-difference element tangent", testElementTangentByFiniteDifference},
        {"element versus dense tangent", testElementResponseAgainstDirectDenseTangent},
    };
    int failures = 0;
    for (const std::pair<std::string, std::function<void()>>& test : tests) {
        try {
            test.second();
            std::cout << "[PASS] " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << test.first << ": " << error.what() << '\n';
        }
    }
    if (failures != 0) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All IGA core tests passed\n";
    return 0;
}
