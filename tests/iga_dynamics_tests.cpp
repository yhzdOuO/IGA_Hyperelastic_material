#include "mechanics/GeneralizedAlphaDynamics.h"
#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <Eigen/Eigenvalues>
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
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> openUniformKnots(int degree, int elements)
{
    std::vector<double> knots;
    for (int i = 0; i <= degree; ++i) knots.push_back(0.0);
    for (int span = 1; span < elements; ++span) {
        knots.push_back(static_cast<double>(span) / elements);
    }
    for (int i = 0; i <= degree; ++i) knots.push_back(1.0);
    return knots;
}

std::vector<double> greville(const std::vector<double>& knots, int degree)
{
    const int count = static_cast<int>(knots.size()) - degree - 1;
    std::vector<double> values(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        for (int j = 1; j <= degree; ++j) values[static_cast<std::size_t>(i)] += knots[i + j];
        values[static_cast<std::size_t>(i)] /= degree;
    }
    return values;
}

TensorProductNurbsSpace makeBox(const Eigen::Vector3d& dimensions,
                                const std::array<int, 3>& degrees,
                                const std::array<int, 3>& elements)
{
    std::array<std::vector<double>, 3> knots{
        openUniformKnots(degrees[0], elements[0]),
        openUniformKnots(degrees[1], elements[1]),
        openUniformKnots(degrees[2], elements[2])};
    std::array<std::vector<double>, 3> coordinates{
        greville(knots[0], degrees[0]), greville(knots[1], degrees[1]),
        greville(knots[2], degrees[2])};
    ControlPointMatrix points(static_cast<int>(coordinates[0].size()
                                               * coordinates[1].size()
                                               * coordinates[2].size()), 3);
    int row = 0;
    for (double z : coordinates[2]) {
        for (double y : coordinates[1]) {
            for (double x : coordinates[0]) {
                points.row(row++) << dimensions(0) * x, dimensions(1) * y,
                                     dimensions(2) * z;
            }
        }
    }
    const int total = static_cast<int>(points.rows());
    return TensorProductNurbsSpace(
        {KnotVector(std::move(knots[0])), KnotVector(std::move(knots[1])),
         KnotVector(std::move(knots[2]))}, degrees, std::move(points),
        Eigen::VectorXd::Ones(total));
}

iga::TimeDependentLoad zeroLoad(int dofs)
{
    return [dofs](double) { return Eigen::VectorXd::Zero(dofs); };
}

void testParametersAndMass()
{
    const iga::GeneralizedAlphaParameters midpoint =
        iga::generalizedAlphaFromSpectralRadius(1.0);
    require(std::abs(midpoint.alphaM - 0.5) < 1e-15
                && std::abs(midpoint.alphaF - 0.5) < 1e-15
                && std::abs(midpoint.gamma - 0.5) < 1e-15
                && std::abs(midpoint.beta - 0.25) < 1e-15,
            "rho_infinity=1 parameter table mismatch");
    const iga::GeneralizedAlphaParameters standard =
        iga::generalizedAlphaFromSpectralRadius(0.8);
    require(std::abs(standard.alphaM - 2.0 / 3.0) < 1e-15
                && std::abs(standard.alphaF - 5.0 / 9.0) < 1e-15
                && std::abs(standard.gamma - 11.0 / 18.0) < 1e-15
                && std::abs(standard.beta - 25.0 / 81.0) < 1e-15,
            "rho_infinity=0.8 parameter table mismatch");

    const TensorProductNurbsSpace space = makeBox({10.0, 1.0, 1.0},
                                                   {3, 2, 2}, {4, 1, 1});
    const iga::SparseMatrix mass = iga::assembleConsistentMass(space, 1.0);
    const iga::MassDiagnostics diagnostics = iga::diagnoseMassMatrix(space, 1.0, mass);
    require(std::abs(diagnostics.integratedMass - 10.0) < 2e-12,
            "consistent total mass must equal 10");
    require(diagnostics.relativeSymmetryError < 1e-14,
            "consistent mass matrix must be symmetric");
    require(diagnostics.positiveDefinite, "consistent mass matrix must be SPD");
}

void testEffectiveTangentFiniteDifference()
{
    const TensorProductNurbsSpace space = makeBox({2.0, 1.0, 1.0},
                                                   {1, 1, 1}, {2, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(4.0, 10.0);
    const iga::SparseMatrix mass = iga::assembleConsistentMass(space, 1.3);
    const int dofs = iga::numDofs(space);
    iga::DynamicState previous;
    previous.time = 0.17;
    previous.displacements.resize(dofs);
    previous.velocities.resize(dofs);
    previous.accelerations.resize(dofs);
    Eigen::VectorXd nextAcceleration(dofs);
    for (int i = 0; i < dofs; ++i) {
        previous.displacements(i) = 0.002 * std::sin(0.31 * (i + 1));
        previous.velocities(i) = 0.01 * std::cos(0.27 * (i + 1));
        previous.accelerations(i) = 0.03 * std::sin(0.19 * (i + 1));
        nextAcceleration(i) = 0.025 * std::cos(0.23 * (i + 1));
    }
    const iga::GeneralizedAlphaParameters parameters =
        iga::generalizedAlphaFromSpectralRadius(0.5);
    const double dt = 0.35;
    const iga::DynamicResidualResult analytic = iga::evaluateDynamicResidual(
        space, material, mass, previous, nextAcceleration, dt, parameters,
        zeroLoad(dofs));
    const Eigen::MatrixXd tangent(analytic.effectiveTangent);
    Eigen::MatrixXd finiteDifference(dofs, dofs);
    const double epsilon = 2.0e-6;
    for (int column = 0; column < dofs; ++column) {
        Eigen::VectorXd plus = nextAcceleration;
        Eigen::VectorXd minus = nextAcceleration;
        plus(column) += epsilon;
        minus(column) -= epsilon;
        finiteDifference.col(column) =
            (iga::evaluateDynamicResidual(space, material, mass, previous, plus,
                                          dt, parameters, zeroLoad(dofs)).residual
             - iga::evaluateDynamicResidual(space, material, mass, previous, minus,
                                             dt, parameters, zeroLoad(dofs)).residual)
            / (2.0 * epsilon);
    }
    const double error = (finiteDifference - tangent).norm() / tangent.norm();
    require(error < 2e-9, "effective tangent finite-difference error="
                              + std::to_string(error));
}

void testZeroStateAndFreeTranslation()
{
    const TensorProductNurbsSpace space = makeBox({2.0, 1.0, 1.0},
                                                   {2, 2, 2}, {2, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(10.0, 20.0);
    const int dofs = iga::numDofs(space);
    iga::DynamicsSettings settings;
    settings.timeStep = 0.025;
    settings.numberOfSteps = 20;
    settings.absoluteTolerance = 1e-11;
    settings.printIterations = false;
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(dofs);
    const iga::DynamicsResult stationary = iga::solveDynamics(
        space, material, 1.0, zero, zero, {}, zeroLoad(dofs), settings);
    require(stationary.states.back().displacements.isZero(0.0),
            "zero state must remain exactly zero");

    Eigen::VectorXd velocity(dofs);
    const Eigen::Vector3d rigidVelocity(0.3, -0.2, 0.1);
    for (int node = 0; node < dofs / 3; ++node) {
        velocity.segment<3>(3 * node) = rigidVelocity;
    }
    const iga::DynamicsResult translating = iga::solveDynamics(
        space, material, 1.0, zero, velocity, {}, zeroLoad(dofs), settings);
    const iga::DynamicState& final = translating.states.back();
    double maxError = 0.0;
    for (int node = 0; node < dofs / 3; ++node) {
        maxError = std::max(maxError,
            (final.displacements.segment<3>(3 * node)
             - final.time * rigidVelocity).cwiseAbs().maxCoeff());
    }
    require(maxError < 2e-12, "free-translation error=" + std::to_string(maxError));
    require(final.accelerations.cwiseAbs().maxCoeff() < 2e-11,
            "free translation must have zero acceleration");
    require(static_cast<int>(translating.history.size()) == settings.numberOfSteps,
            "exact predictor should converge on the first iteration of every step");
}

void testInitialAccelerationAndEnergy()
{
    const TensorProductNurbsSpace space = makeBox({2.0, 1.0, 1.0},
                                                   {2, 2, 2}, {2, 1, 1});
    const iga::SaintVenantKirchhoffMaterial material(0.0, 1000.0);
    const int dofs = iga::numDofs(space);
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    const Eigen::VectorXd zero = Eigen::VectorXd::Zero(dofs);
    const Eigen::VectorXd body = iga::assembleBodyForce(
        space, 2.0, Eigen::Vector3d(0.0, 0.0, -3.0));
    const iga::SparseMatrix mass = iga::assembleConsistentMass(space, 2.0);
    const iga::TimeDependentLoad load = [body](double) { return body; };
    const Eigen::VectorXd acceleration = iga::computeInitialAcceleration(
        space, material, mass, zero, fixed, load);
    const Eigen::VectorXd equilibrium = mass * acceleration - body;
    require(iga::freeDofNorm(equilibrium, fixed) < 2e-12,
            "initial acceleration must satisfy the t=0 equation");
    for (int dof : fixed) require(acceleration(dof) == 0.0,
                                  "constrained initial acceleration must be zero");
    require(iga::computeStrainEnergy(space, zero, material) == 0.0,
            "undeformed state must have zero strain energy");
}

Eigen::VectorXd axialModeDisplacement(const TensorProductNurbsSpace& space,
                                      double amplitude)
{
    Eigen::VectorXd displacement = Eigen::VectorXd::Zero(iga::numDofs(space));
    for (int node = 0; node < space.controlPoints().rows(); ++node) {
        const double x = space.controlPoints()(node, 0);
        displacement(3 * node) = amplitude * std::sin(0.5 * M_PI * x / 10.0);
    }
    return displacement;
}

iga::DynamicsResult solveAxialVibration(const TensorProductNurbsSpace& space,
                                        const Eigen::VectorXd& initialDisplacement,
                                        double timeStep, int steps)
{
    const double youngsModulus = 1.0e7;
    const iga::SaintVenantKirchhoffMaterial material(0.0, 0.5 * youngsModulus);
    iga::DynamicsSettings settings;
    settings.timeStep = timeStep;
    settings.numberOfSteps = steps;
    settings.spectralRadius = 1.0;
    settings.relativeTolerance = 1.0e-10;
    settings.absoluteTolerance = 1.0e-13;
    settings.printIterations = false;
    return iga::solveDynamics(
        space, material, 1.0, initialDisplacement,
        Eigen::VectorXd::Zero(iga::numDofs(space)), iga::collectXiMinFaceDofs(space),
        zeroLoad(iga::numDofs(space)), settings);
}

void testAxialFrequencyAndSecondOrderAccuracy()
{
    const TensorProductNurbsSpace space = makeBox({10.0, 1.0, 1.0},
                                                   {3, 1, 1}, {8, 1, 1});
    const double exactPeriod = 4.0 * 10.0 / std::sqrt(1.0e7);
    const Eigen::VectorXd approximateMode = axialModeDisplacement(space, 1.0e-6);
    const double fineDt = exactPeriod / 80.0;
    const iga::DynamicsResult frequency = solveAxialVibration(
        space, approximateMode, fineDt, 200);
    std::vector<double> fallingCrossings;
    double previousTip = iga::evaluateDisplacement(
        space, frequency.states.front().displacements, {1.0, 0.5, 0.5})(0);
    for (std::size_t i = 1; i < frequency.states.size(); ++i) {
        const double tip = iga::evaluateDisplacement(
            space, frequency.states[i].displacements, {1.0, 0.5, 0.5})(0);
        if (previousTip >= 0.0 && tip < 0.0) {
            const double fraction = previousTip / (previousTip - tip);
            fallingCrossings.push_back(frequency.states[i - 1].time + fraction * fineDt);
        }
        previousTip = tip;
    }
    require(fallingCrossings.size() >= 2, "axial vibration did not complete one period");
    const double numericalPeriod = fallingCrossings[1] - fallingCrossings[0];
    require(std::abs(numericalPeriod - exactPeriod) / exactPeriod < 0.01,
            "axial vibration period error="
                + std::to_string(std::abs(numericalPeriod - exactPeriod) / exactPeriod));
    const double initialEnergy = frequency.states.front().kineticEnergy
        + frequency.states.front().strainEnergy;
    double maximumEnergyDrift = 0.0;
    for (const iga::DynamicState& state : frequency.states) {
        maximumEnergyDrift = std::max(maximumEnergyDrift,
            std::abs(state.kineticEnergy + state.strainEnergy - initialEnergy)
                / initialEnergy);
    }
    require(maximumEnergyDrift < 2e-8,
            "midpoint energy drift in the linear regime="
                + std::to_string(maximumEnergyDrift));

    // Build the discrete axial eigenmode so the linear semi-discrete solution is known
    // exactly; this isolates temporal error from spatial modal contamination.
    const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    const iga::SparseMatrix massSparse = iga::assembleConsistentMass(space, 1.0);
    const iga::SparseMatrix stiffnessSparse = iga::assembleInternalSystem(
        space, Eigen::VectorXd::Zero(iga::numDofs(space)), material).tangentStiffness;
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    std::vector<bool> isFixed(static_cast<std::size_t>(iga::numDofs(space)), false);
    for (int dof : fixed) isFixed[static_cast<std::size_t>(dof)] = true;
    std::vector<int> free;
    for (int dof = 0; dof < iga::numDofs(space); ++dof) {
        if (!isFixed[static_cast<std::size_t>(dof)]) free.push_back(dof);
    }
    Eigen::MatrixXd mass(free.size(), free.size());
    Eigen::MatrixXd stiffness(free.size(), free.size());
    const Eigen::MatrixXd fullMass(massSparse);
    const Eigen::MatrixXd fullStiffness(stiffnessSparse);
    for (std::size_t i = 0; i < free.size(); ++i) {
        for (std::size_t j = 0; j < free.size(); ++j) {
            mass(i, j) = fullMass(free[i], free[j]);
            stiffness(i, j) = fullStiffness(free[i], free[j]);
        }
    }
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> modes(stiffness, mass);
    require(modes.info() == Eigen::Success, "discrete modal analysis failed");
    const double targetEigenvalue = std::pow(2.0 * M_PI / exactPeriod, 2);
    int axialMode = 0;
    (modes.eigenvalues().array() - targetEigenvalue).abs().minCoeff(&axialMode);
    const double omega = std::sqrt(modes.eigenvalues()(axialMode));
    Eigen::VectorXd discreteMode = Eigen::VectorXd::Zero(iga::numDofs(space));
    for (std::size_t i = 0; i < free.size(); ++i) {
        discreteMode(free[i]) = modes.eigenvectors()(static_cast<int>(i), axialMode);
    }
    discreteMode *= 1.0e-6 / discreteMode.cwiseAbs().maxCoeff();

    const double discretePeriod = 2.0 * M_PI / omega;
    const double finalTime = 0.37 * discretePeriod;
    const int coarseSteps = 12;
    const iga::DynamicsResult coarse = solveAxialVibration(
        space, discreteMode, finalTime / coarseSteps, coarseSteps);
    const iga::DynamicsResult medium = solveAxialVibration(
        space, discreteMode, finalTime / (2 * coarseSteps), 2 * coarseSteps);
    const iga::DynamicsResult fine = solveAxialVibration(
        space, discreteMode, finalTime / (4 * coarseSteps), 4 * coarseSteps);
    const Eigen::VectorXd exact = std::cos(omega * finalTime) * discreteMode;
    const double coarseError = (coarse.states.back().displacements - exact).norm();
    const double mediumError = (medium.states.back().displacements - exact).norm();
    const double fineError = (fine.states.back().displacements - exact).norm();
    require(coarseError / mediumError > 3.5 && mediumError / fineError > 3.5,
            "time refinement must exhibit second-order convergence; ratio="
                + std::to_string(coarseError / mediumError) + ", "
                + std::to_string(mediumError / fineError));
}

void testCantileverStepLoadAmplification()
{
    const TensorProductNurbsSpace space = makeBox({10.0, 1.0, 1.0},
                                                   {3, 3, 3}, {16, 2, 2});
    const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    const std::vector<int> fixed = iga::collectXiMinFaceDofs(space);
    const Eigen::VectorXd force = iga::assembleBodyForce(
        space, 1.0, Eigen::Vector3d(0.0, 0.0, -1.0));
    iga::NewtonSettings staticSettings;
    staticSettings.printIterations = false;
    staticSettings.relativeTolerance = 1e-11;
    const iga::NewtonResult staticResult = iga::solveDeadLoad(
        space, material, force, fixed, staticSettings);
    const double staticTip = -iga::evaluateDisplacement(
        space, staticResult.displacements, {1.0, 0.5, 0.5})(2);

    iga::DynamicsSettings settings;
    settings.timeStep = 4.9e-3;
    // Run well beyond the third bending period; this guards against the delayed
    // instability that previously appeared around steps 125-135.
    settings.numberOfSteps = 180;
    settings.spectralRadius = 0.8;
    settings.printIterations = false;
    const iga::TimeDependentLoad load = [force](double) { return force; };
    const iga::DynamicsResult dynamics = iga::solveDynamics(
        space, material, 1.0, Eigen::VectorXd::Zero(iga::numDofs(space)),
        Eigen::VectorXd::Zero(iga::numDofs(space)), fixed, load, settings);
    double maximumTip = 0.0;
    for (const iga::DynamicState& state : dynamics.states) {
        const double tip = -iga::evaluateDisplacement(
            space, state.displacements, {1.0, 0.5, 0.5})(2);
        require(std::isfinite(tip) && std::abs(tip) < 1.0e-2,
                "cantilever response became unbounded at t="
                    + std::to_string(state.time));
        require(std::isfinite(state.kineticEnergy)
                    && std::isfinite(state.strainEnergy),
                "cantilever energy became non-finite at t="
                    + std::to_string(state.time));
        maximumTip = std::max(maximumTip, tip);
    }
    const double amplification = maximumTip / staticTip;
    require(amplification > 1.8 && amplification < 2.15,
            "cantilever step-load amplification must be about two; value="
                + std::to_string(amplification));
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"generalized-alpha parameters and mass", testParametersAndMass},
        {"effective tangent finite difference", testEffectiveTangentFiniteDifference},
        {"zero state and free translation", testZeroStateAndFreeTranslation},
        {"initial acceleration and energy", testInitialAccelerationAndEnergy},
        {"axial frequency and second-order accuracy",
         testAxialFrequencyAndSecondOrderAccuracy},
        {"cantilever step-load amplification", testCantileverStepLoadAmplification},
    };
    int failures = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "[PASS] " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << test.first << ": " << error.what() << '\n';
        }
    }
    if (failures != 0) return 1;
    std::cout << "All dynamics tests passed\n";
    return 0;
}
