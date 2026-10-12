#include "mechanics/GlobalNonlinearSystem.h"
#include "mechanics/NeumannBoundaryLoad.h"
#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <Eigen/SparseCholesky>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
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

std::vector<double> openUniformKnots(int degree, int elements)
{
    std::vector<double> knots;
    knots.reserve(static_cast<std::size_t>(2 * (degree + 1) + elements - 1));
    for (int i = 0; i <= degree; ++i) {
        knots.push_back(0.0);
    }
    for (int span = 1; span < elements; ++span) {
        knots.push_back(static_cast<double>(span) / elements);
    }
    for (int i = 0; i <= degree; ++i) {
        knots.push_back(1.0);
    }
    return knots;
}

std::vector<double> grevilleAbscissae(const std::vector<double>& knots,
                                      int degree)
{
    const int count = static_cast<int>(knots.size()) - degree - 1;
    std::vector<double> result(static_cast<std::size_t>(count), 0.0);
    for (int i = 0; i < count; ++i) {
        double sum = 0.0;
        for (int j = 1; j <= degree; ++j) {
            sum += knots[static_cast<std::size_t>(i + j)];
        }
        result[static_cast<std::size_t>(i)] = sum / degree;
    }
    return result;
}

TensorProductNurbsSpace makeBox(
    const Eigen::Vector3d& dimensions,
    const std::array<int, 3>& degrees,
    const std::array<int, 3>& elements)
{
    std::array<std::vector<double>, 3> knotValues{
        openUniformKnots(degrees[0], elements[0]),
        openUniformKnots(degrees[1], elements[1]),
        openUniformKnots(degrees[2], elements[2])};
    std::array<std::vector<double>, 3> greville{
        grevilleAbscissae(knotValues[0], degrees[0]),
        grevilleAbscissae(knotValues[1], degrees[1]),
        grevilleAbscissae(knotValues[2], degrees[2])};
    const int n0 = static_cast<int>(greville[0].size());
    const int n1 = static_cast<int>(greville[1].size());
    const int n2 = static_cast<int>(greville[2].size());
    ControlPointMatrix points(n0 * n1 * n2, 3);
    int global = 0;
    for (int i2 = 0; i2 < n2; ++i2) {
        for (int i1 = 0; i1 < n1; ++i1) {
            for (int i0 = 0; i0 < n0; ++i0) {
                points.row(global++) << dimensions(0) * greville[0][static_cast<std::size_t>(i0)],
                                        dimensions(1) * greville[1][static_cast<std::size_t>(i1)],
                                        dimensions(2) * greville[2][static_cast<std::size_t>(i2)];
            }
        }
    }
    const int total = static_cast<int>(points.rows());
    return TensorProductNurbsSpace(
        {KnotVector(std::move(knotValues[0])),
         KnotVector(std::move(knotValues[1])),
         KnotVector(std::move(knotValues[2]))},
        degrees, std::move(points), Eigen::VectorXd::Ones(total));
}

Eigen::Vector3d sumByComponent(const Eigen::VectorXd& vector)
{
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    for (int basis = 0; basis < vector.size() / 3; ++basis) {
        result += vector.segment<3>(3 * basis);
    }
    return result;
}

void testReferenceFaceTractionLoads()
{
    const Eigen::Vector3d dimensions(2.0, 3.0, 4.0);
    const TensorProductNurbsSpace space =
        makeBox(dimensions, {2, 3, 1}, {2, 2, 2});
    const std::array<double, 3> faceAreas{
        dimensions(1) * dimensions(2),
        dimensions(0) * dimensions(2),
        dimensions(0) * dimensions(1)};
    const std::array<int, 3> counts{
        space.numControlPoints(0),
        space.numControlPoints(1),
        space.numControlPoints(2)};

    // ---------- 功能块 1：六个参数面都应满足离散合力等于牵引力乘参考面积 ----------
    // 使用不同次数和不同 span 数，避免测试只覆盖三向完全相同的特殊网格。
    for (int direction = 0; direction < 3; ++direction) {
        for (int side = 0; side < 2; ++side) {
            int callbackCalls = 0;
            const iga::TimeDependentLoad load = iga::makeReferenceFaceTractionLoad(
                space, direction, side,
                [&callbackCalls](double time) {
                    ++callbackCalls;
                    return Eigen::Vector3d(1.0 + time, -2.0 * time, 0.5);
                });
            require(callbackCalls == 0,
                    "traction callback must not run while spatial weights are preintegrated");
            const double time = 0.25;
            const Eigen::Vector3d traction(1.25, -0.5, 0.5);
            const Eigen::VectorXd force = load(time);
            require(callbackCalls == 1,
                    "each global load evaluation must call the traction callback once");
            const Eigen::Vector3d expected = faceAreas[static_cast<std::size_t>(direction)]
                * traction;
            require((sumByComponent(force) - expected).norm()
                        < 2.0e-12 * expected.norm(),
                    "reference face resultant must equal traction times reference area");

            // 功能块 1.1：开放 knot vector 的端面载荷只能分配到相应边界控制层。
            for (int basis = 0; basis < space.controlPoints().rows(); ++basis) {
                const std::array<int, 3> indices{
                    basis % counts[0],
                    (basis / counts[0]) % counts[1],
                    basis / (counts[0] * counts[1])};
                const int boundaryIndex = side == 0
                    ? 0 : counts[static_cast<std::size_t>(direction)] - 1;
                if (indices[static_cast<std::size_t>(direction)] != boundaryIndex) {
                    require(force.segment<3>(3 * basis).norm() < 1.0e-14,
                            "face traction leaked to a control point outside the selected side");
                }
            }
        }
    }

    // ---------- 功能块 2：在接口边界拒绝无效面编号和无效回调 ----------
    bool rejectedDirection = false;
    try {
        iga::makeReferenceFaceTractionLoad(
            space, 3, 0, [](double) { return Eigen::Vector3d::Zero(); });
    } catch (const std::invalid_argument&) {
        rejectedDirection = true;
    }
    require(rejectedDirection, "invalid face direction must be rejected");

    bool rejectedSide = false;
    try {
        iga::makeReferenceFaceTractionLoad(
            space, 0, -1, [](double) { return Eigen::Vector3d::Zero(); });
    } catch (const std::invalid_argument&) {
        rejectedSide = true;
    }
    require(rejectedSide, "invalid face side must be rejected");

    bool rejectedEmptyCallback = false;
    try {
        iga::makeReferenceFaceTractionLoad(
            space, 0, 0, iga::ReferenceTraction{});
    } catch (const std::invalid_argument&) {
        rejectedEmptyCallback = true;
    }
    require(rejectedEmptyCallback, "empty traction callback must be rejected");
}

void testBodyForceAndDofOrdering()
{
    const TensorProductNurbsSpace space =
        makeBox({2.0, 3.0, 4.0}, {2, 2, 2}, {2, 1, 1});
    const Eigen::Vector3d acceleration(1.5, -2.0, 0.25);
    const double density = 7.0;
    const Eigen::VectorXd force =
        iga::assembleBodyForce(space, density, acceleration);
    const Eigen::Vector3d expected = density * (2.0 * 3.0 * 4.0) * acceleration;
    const double forceBalanceError = (sumByComponent(force) - expected).norm();
    require(forceBalanceError / expected.norm() < 1e-12,
            "consistent body force must equal density*volume*acceleration; error="
                + std::to_string(forceBalanceError));
    require(iga::numDofs(space) == 3 * space.controlPoints().rows(),
            "global DOF count must be 3*Ncp");

    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    require(static_cast<int>(fixed.size())
                == 3 * space.numControlPoints(1) * space.numControlPoints(2),
            "xi-min face constraint count is incorrect");
    for (std::size_t index = 0; index < fixed.size(); index += 3) {
        require(fixed[index + 1] == fixed[index] + 1
                    && fixed[index + 2] == fixed[index] + 2,
                "constraint DOFs must use interleaved 3*A+i ordering");
    }
}

void testAssemblyBoundaryConditionsAndGlobalFiniteDifference()
{
    const TensorProductNurbsSpace space =
        makeBox({2.0, 1.0, 1.0}, {1, 1, 1}, {2, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(4.0, 2.0);
    Eigen::VectorXd displacement(iga::numDofs(space));
    for (int dof = 0; dof < displacement.size(); ++dof) {
        displacement(dof) = 0.004 * std::sin(0.31 * (dof + 1));
    }
    const iga::GlobalAssemblyResult assembly =
        iga::assembleInternalSystem(space, displacement, material);
    const Eigen::MatrixXd dense = Eigen::MatrixXd(assembly.tangentStiffness);
    require((dense - dense.transpose()).norm() / dense.norm() < 1e-13,
            "assembled tangent must be symmetric");

    const double step = 2.0e-6;
    Eigen::MatrixXd finiteDifference(dense.rows(), dense.cols());
    for (int column = 0; column < dense.cols(); ++column) {
        Eigen::VectorXd plus = displacement;
        Eigen::VectorXd minus = displacement;
        plus(column) += step;
        minus(column) -= step;
        finiteDifference.col(column) =
            (iga::assembleInternalSystem(space, plus, material).internalForce
             - iga::assembleInternalSystem(space, minus, material).internalForce)
            / (2.0 * step);
    }
    require((finiteDifference - dense).norm() / dense.norm() < 2e-10,
            "global finite-difference tangent check failed");

    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    const iga::ConstrainedLinearSystem constrained =
        iga::applyHomogeneousDirichlet(assembly.tangentStiffness,
                                      -assembly.internalForce, fixed);
    const Eigen::MatrixXd constrainedDense = Eigen::MatrixXd(constrained.matrix);
    require((constrainedDense - constrainedDense.transpose()).norm() < 1e-13,
            "Dirichlet processing must preserve symmetry");
    Eigen::SimplicialLDLT<iga::SparseMatrix> solver;
    solver.compute(constrained.matrix);
    require(solver.info() == Eigen::Success,
            "fully fixed xi face must remove all rigid-body modes");
}

iga::NewtonSettings quietNewtonSettings(int loadSteps = 1)
{
    iga::NewtonSettings settings;
    settings.loadSteps = loadSteps;
    settings.maxIterations = 15;
    settings.relativeTolerance = 1e-10;
    settings.absoluteTolerance = 1e-13;
    settings.printIterations = false;
    return settings;
}

void testZeroLoadAndAxialRod()
{
    const TensorProductNurbsSpace space =
        makeBox({10.0, 1.0, 1.0}, {3, 3, 3}, {8, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(iga::numDofs(space));
    const iga::NewtonResult zeroResult =
        iga::solveDeadLoad(space, material, zero, fixed, quietNewtonSettings());
    require(zeroResult.displacements.isZero(0.0)
                && zeroResult.history.size() == 1,
            "zero load must converge immediately to zero displacement");

    const Eigen::VectorXd force =
        iga::assembleBodyForce(space, 1.0, Eigen::Vector3d(1.0, 0.0, 0.0));
    const iga::NewtonResult result =
        iga::solveDeadLoad(space, material, force, fixed, quietNewtonSettings());
    const Eigen::Vector3d tip = iga::evaluateDisplacement(
        space, result.displacements, {1.0, 0.5, 0.5});
    const double exactTip = 5.0e-6;
    require(std::abs(tip(0) - exactTip) / exactTip < 1e-5,
            "small-load axial rod must match the exact quadratic solution");
    require(std::max(std::abs(tip(1)), std::abs(tip(2))) < 1e-11,
            "nu=0 axial rod must have zero transverse displacement");

    const Eigen::Vector3d reaction = sumByComponent(result.residual);
    require((reaction + sumByComponent(force)).norm() < 2e-8,
            "support reactions must balance the total body force");
}

void testNonlinearAxialRod()
{
    const TensorProductNurbsSpace space =
        makeBox({10.0, 1.0, 1.0}, {3, 3, 3}, {8, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    const Eigen::VectorXd force =
        iga::assembleBodyForce(space, 1.0, Eigen::Vector3d(1.0e4, 0.0, 0.0));
    const iga::NewtonResult result =
        iga::solveDeadLoad(space, material, force, fixed, quietNewtonSettings());
    const double tip = iga::evaluateDisplacement(
        space, result.displacements, {1.0, 0.5, 0.5})(0);
    const double exactTip = 4.950974524113e-2;
    require(std::abs(tip - exactTip) / exactTip < 2e-5,
            "nonlinear axial rod must match the SVK reference solution");
    require(std::abs(tip - 5.0e-2) > 4.0e-4,
            "nonlinear axial result must be distinguishable from the linear prediction");
}

void testCantileverLinearAndNonlinearBehavior()
{
    const TensorProductNurbsSpace space =
        makeBox({10.0, 1.0, 1.0}, {3, 3, 3}, {16, 2, 2});
    const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);

    const Eigen::VectorXd smallForce =
        iga::assembleBodyForce(space, 1.0, Eigen::Vector3d(0.0, 0.0, -1.0));
    const iga::NewtonResult small =
        iga::solveDeadLoad(space, material, smallForce, fixed, quietNewtonSettings());
    const double smallTip = -iga::evaluateDisplacement(
        space, small.displacements, {1.0, 0.5, 0.5})(2);
    require(std::abs(smallTip - 1.5e-3) / 1.5e-3 < 0.05,
            "linear cantilever tip displacement must be within 5% of Euler-Bernoulli");

    const Eigen::VectorXd largeForce = 500.0 * smallForce;
    const iga::NewtonResult large =
        iga::solveDeadLoad(space, material, largeForce, fixed, quietNewtonSettings(1));
    const double largeTip = -iga::evaluateDisplacement(
        space, large.displacements, {1.0, 0.5, 0.5})(2);
    std::cout << "  cantilever small tip=" << smallTip
              << ", large tip=" << largeTip << ", residuals:";
    for (const iga::NewtonIterationRecord& record : large.history) {
        std::cout << ' ' << record.relativeResidual;
    }
    std::cout << '\n';
    // Compare with this same 3-D discretization's linear extrapolation. Comparing
    // against 500*the Euler-Bernoulli value mixes geometric stiffening with the
    // mesh's small shear/solid-model offset and can reject a correctly stiffened result.
    const double linearExtrapolation = 500.0 * smallTip;
    require(largeTip > 0.0 && largeTip < 0.999 * linearExtrapolation,
            "large-load cantilever must be stiffer than its discrete linear extrapolation; tip="
                + std::to_string(largeTip));

    std::vector<double> lastStepResiduals;
    for (const iga::NewtonIterationRecord& record : large.history) {
        if (record.loadStep == 1) {
            lastStepResiduals.push_back(record.relativeResidual);
        }
    }
    require(lastStepResiduals.size() >= 4 && lastStepResiduals.size() <= 7,
            "nonlinear cantilever should require roughly 4-5 Newton iterations");
    const std::size_t count = lastStepResiduals.size();
    require(lastStepResiduals[count - 1] < 1e-9
                && lastStepResiduals[count - 2] < 1e-4,
            "nonlinear cantilever must display the final quadratic residual drop");
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"reference face traction loads", testReferenceFaceTractionLoads},
        {"body force and DOF ordering", testBodyForceAndDofOrdering},
        {"assembly, BC, and global finite difference",
         testAssemblyBoundaryConditionsAndGlobalFiniteDifference},
        {"zero load and exact axial rod", testZeroLoadAndAxialRod},
        {"nonlinear axial rod", testNonlinearAxialRod},
        {"linear and nonlinear cantilever", testCantileverLinearAndNonlinearBehavior},
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
        std::cerr << failures << " global test(s) failed\n";
        return 1;
    }
    std::cout << "All stage-two global tests passed\n";
    return 0;
}
