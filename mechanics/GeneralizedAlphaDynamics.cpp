#include "mechanics/GeneralizedAlphaDynamics.h"

#include <Eigen/SparseCholesky>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace iga {
namespace {

// 动力学入口会反复接收用户初值和载荷回调结果；统一在边界处检查
// N_dof 维度与有限性，避免 NaN 进入稀疏分解后才以 NumericalIssue 暴露。
void validateStateVector(const Eigen::VectorXd& vector, int size, const char* name)
{
    if (vector.size() != size || !vector.allFinite()) {
        throw std::invalid_argument(std::string(name) + " has invalid size or values");
    }
}

void requireHomogeneousConstraints(const Eigen::VectorXd& vector,
                                   const std::vector<int>& constrainedDofs,
                                   const char* name)
{
    // 本阶段只支持齐次位移约束，因此 d、v、a 在受约束位置都必须为零。
    for (int dof : constrainedDofs) {
        if (dof < 0 || dof >= vector.size()) {
            throw std::out_of_range("constrained DOF is outside the global vector");
        }
        if (std::abs(vector(dof)) > 1.0e-14) {
            throw std::invalid_argument(std::string(name)
                                        + " must satisfy homogeneous constraints");
        }
    }
}

Eigen::VectorXd checkedLoad(const TimeDependentLoad& load, double time, int size)
{
    // 空 std::function 是合法的零载荷表达，便于自由振动和刚体平动测试。
    if (!load) {
        return Eigen::VectorXd::Zero(size);
    }
    Eigen::VectorXd value = load(time);
    validateStateVector(value, size, "external load");
    return value;
}

DynamicState makeState(double time,
                       Eigen::VectorXd displacements,
                       Eigen::VectorXd velocities,
                       Eigen::VectorXd accelerations,
                       const SparseMatrix& massMatrix,
                       const nurbs::TensorProductNurbsSpace& space,
                       const HyperelasticMaterial& material)
{
    // 只为已经收敛的物理时刻构造状态；中间构形不进入历史，也不交给渲染。
    DynamicState state;
    state.time = time;
    state.displacements = std::move(displacements);
    state.velocities = std::move(velocities);
    state.accelerations = std::move(accelerations);
    state.kineticEnergy = 0.5 * state.velocities.dot(massMatrix * state.velocities);
    state.strainEnergy = computeStrainEnergy(space, state.displacements, material);
    return state; //当前时间步收敛状态的一份记录
}

} // namespace

GeneralizedAlphaParameters generalizedAlphaFromSpectralRadius(double spectralRadius)
{
    if (!std::isfinite(spectralRadius) || spectralRadius < 0.0
        || spectralRadius > 1.0) {
        throw std::invalid_argument("spectral radius must lie in [0, 1]");
    }
    // 本项目采用 alpha=1 表示全隐式的约定。以下四式必须作为一组使用，
    // 否则即使求解稳定，时间精度也会从二阶退化为一阶。
    GeneralizedAlphaParameters result;
    result.spectralRadius = spectralRadius;
    result.alphaM = (2.0 - spectralRadius) / (1.0 + spectralRadius);
    result.alphaF = 1.0 / (1.0 + spectralRadius);
    result.gamma = 0.5 - result.alphaF + result.alphaM;
    const double secondOrderFactor = 1.0 - result.alphaF + result.alphaM;
    result.beta = 0.25 * secondOrderFactor * secondOrderFactor;
    return result;
}

SparseMatrix assembleConsistentMass(const nurbs::TensorProductNurbsSpace& space,
                                    double referenceDensity)
{
    if (!std::isfinite(referenceDensity) || referenceDensity <= 0.0) {
        throw std::invalid_argument("reference density must be finite and positive");
    }
    // ---------- 功能块 1：逐单元积分标量一致质量 Me ----------
    // Me_ab=integral(rho0*R_a*R_b)dOmega0，大小为 n_en x n_en。
    const int nodes = space.numElementNodes();
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(space.numElements() * nodes * nodes * 3));
    for (int element = 0; element < space.numElements(); ++element) {
        // &逐element处理
        Eigen::MatrixXd scalarMass = Eigen::MatrixXd::Zero(nodes, nodes);
        for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
            for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
                for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                    const ReferencePointData reference = evaluateReferencePoint(
                        space, element, {q0, q1, q2});
                    // rankUpdate 只写下三角；权重已包含参考几何和母单元 Jacobian。
                    scalarMass.selfadjointView<Eigen::Lower>().rankUpdate(
                        reference.shapeValues,
                        referenceDensity * reference.integrationWeight); // &计算局部的质量矩阵M+=α*R*R^T
                }
            }
        }
        // 必须显式镜像下三角，否则随后散射会得到非对称/非正定的全局 M。
        scalarMass = scalarMass.selfadjointView<Eigen::Lower>();

        // ---------- 功能块 2：Me tensor I3 散射到全局自由度 ----------
        // 每个标量 Me_ab 只进入 (A,x/y/z)-(B,x/y/z) 三个同分量位置；
        // 不同位移分量之间的质量耦合为零。
        const std::vector<int> connectivity = space.elementConnectivity(element); // &IEN(e, a)=A => IEN_e(a)=A其中IEN_e就是connectivity
        for (int a = 0; a < nodes; ++a) {
            for (int b = 0; b < nodes; ++b) {
                const double value = scalarMass(a, b);
                // &散射到I3并装配
                for (int component = 0; component < 3; ++component) {
                    triplets.emplace_back(3 * connectivity[static_cast<std::size_t>(a)] + component,
                                          3 * connectivity[static_cast<std::size_t>(b)] + component,
                                          value); // &会有重复triplet
                }
            }
        }
    }
    // 重复 triplet 由 Eigen 自动累加，得到 N_dof x N_dof 的完整一致质量矩阵。
    SparseMatrix mass(numDofs(space), numDofs(space));
    mass.setFromTriplets(triplets.begin(), triplets.end());
    mass.makeCompressed();
    return mass; //&返回全局一致质量矩阵
}

MassDiagnostics diagnoseMassMatrix(const nurbs::TensorProductNurbsSpace& space,
                                   double referenceDensity,
                                   const SparseMatrix& massMatrix)
{
    if (massMatrix.rows() != numDofs(space) || massMatrix.cols() != numDofs(space)) {
        throw std::invalid_argument("mass matrix dimensions are inconsistent");
    }
    // ---------- 三个相互独立的质量矩阵自检 ----------
    MassDiagnostics result;
    const Eigen::VectorXd ones = Eigen::VectorXd::Ones(massMatrix.rows());
    // 三个位移分量各包含一份相同总质量，因此 1^T*M*1 需要除以 3。
    result.integratedMass = ones.dot(massMatrix * ones) / 3.0;
    // 用单位 x 加速度的体力积分独立得到 rho0*V0，避免用同一公式自证。
    result.expectedMass = assembleBodyForce(space, referenceDensity,
                                            Eigen::Vector3d::UnitX()).sum();
    result.relativeMassError = std::abs(result.integratedMass - result.expectedMass)
        / std::max(std::abs(result.expectedMass), std::numeric_limits<double>::min());
    const SparseMatrix difference = massMatrix - SparseMatrix(massMatrix.transpose());
    result.relativeSymmetryError = difference.norm()
        / std::max(massMatrix.norm(), std::numeric_limits<double>::min());
    Eigen::SimplicialLLT<SparseMatrix> cholesky;
    cholesky.compute(massMatrix);
    result.positiveDefinite = cholesky.info() == Eigen::Success;
    return result;
}

DynamicResidualResult evaluateDynamicResidual(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const SparseMatrix& massMatrix,
    const DynamicState& previousState,
    const Eigen::VectorXd& nextAcceleration,
    double timeStep,
    const GeneralizedAlphaParameters& parameters,
    const TimeDependentLoad& externalLoad,
    ElementFormulation formulation)
{
    const int dofs = numDofs(space);
    validateStateVector(previousState.displacements, dofs, "previous displacements");
    validateStateVector(previousState.velocities, dofs, "previous velocities");
    validateStateVector(previousState.accelerations, dofs, "previous accelerations");
    validateStateVector(nextAcceleration, dofs, "next acceleration");
    if (!std::isfinite(timeStep) || timeStep <= 0.0) {
        throw std::invalid_argument("time step must be finite and positive");
    }
    // ---------- 功能块 1：由 Newton 未知量 a_{n+1} 重建 d_{n+1} ----------
    // 这里使用完整 Newmark 关系；d、v 的 Newton 增量系数分别是
    // beta*dt^2 和 gamma*dt，不能额外乘 alpha_f。
    const double dt2 = timeStep * timeStep;
    //& d_{n+1}^{i+1} = d_n + Δt * v_n + Δt^2 * 0.5 *((1-2β) * a_n + 2β * a_{n+1}^{i+1})
    const Eigen::VectorXd nextDisplacement = previousState.displacements
        + timeStep * previousState.velocities
        + dt2 * ((0.5 - parameters.beta) * previousState.accelerations
                 + parameters.beta * nextAcceleration);

    // ---------- 功能块 2：在 generalized-alpha 中间时刻求残量 ----------
    DynamicResidualResult result;
    result.displacementAlphaF = previousState.displacements
        + parameters.alphaF * (nextDisplacement - previousState.displacements); //&d_{n+α_f}
    result.accelerationAlphaM = previousState.accelerations
        + parameters.alphaM * (nextAcceleration - previousState.accelerations); //&a_{n+α_m}
    const GlobalAssemblyResult assembly = assembleInternalSystem(
        space, result.displacementAlphaF, material, formulation); //& d_{n+α_f}时的K和f^{int}
    // 位移、内力和外载荷取 t_{n+alpha_f}；惯性项单独取 a_{n+alpha_m}。
    const double evaluationTime = previousState.time + parameters.alphaF * timeStep;
    //& R(d_{n+α_f}, a_{n+α_m}) = M*a_{n+α_m} + f^{int}(d_{n+α_f}) - f^{ext}(t)
    result.residual = massMatrix * result.accelerationAlphaM
        + assembly.internalForce - checkedLoad(externalLoad, evaluationTime, dofs);
    // 链式法则 dR/da_{n+1}。当前材料无速度型阻尼，所以没有 C 项。
    //& M^* = α_m * M + (α_f * β * Δt^2) * K
    result.effectiveTangent = parameters.alphaM * massMatrix
        + (parameters.alphaF * parameters.beta * dt2) * assembly.tangentStiffness;
    result.effectiveTangent.makeCompressed();
    return result; //& M^* 和 R(d_{n+α_f}, a_{n+α_m})
    //& M^* * Δa_{n+1} = R，其中M^*和R会随着i和n变化
}

Eigen::VectorXd computeInitialAcceleration(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const SparseMatrix& massMatrix,
    const Eigen::VectorXd& initialDisplacements,
    const std::vector<int>& constrainedDofs,
    const TimeDependentLoad& externalLoad,
    ElementFormulation formulation)
{
    const int dofs = numDofs(space);
    validateStateVector(initialDisplacements, dofs, "initial displacements");
    requireHomogeneousConstraints(initialDisplacements, constrainedDofs,
                                  "initial displacements");
    // ---------- 初始平衡：M*a0=f_ext(0)-f_int(d0) ----------
    const Eigen::VectorXd internalForce = assembleInternalSystem(
        space, initialDisplacements, material, formulation).internalForce;
    const ConstrainedLinearSystem constrained = applyHomogeneousDirichlet(
        massMatrix, checkedLoad(externalLoad, 0.0, dofs) - internalForce,
        constrainedDofs);
    Eigen::SimplicialLLT<SparseMatrix> solver;
    // M 在施加齐次约束后仍为 SPD，因此使用 LLT，而不是静力切线所需的 LDLT。
    solver.compute(constrained.matrix);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("initial acceleration mass solve factorization failed");
    }
    Eigen::VectorXd acceleration = solver.solve(constrained.rightHandSide);
    if (solver.info() != Eigen::Success || !acceleration.allFinite()) {
        throw std::runtime_error("initial acceleration mass solve failed");
    }
    return acceleration; //初始状态下的加速度
}

DynamicsResult solveDynamics(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    double referenceDensity,
    const Eigen::VectorXd& initialDisplacements,
    const Eigen::VectorXd& initialVelocities,
    const std::vector<int>& constrainedDofs,
    const TimeDependentLoad& externalLoad,
    const DynamicsSettings& settings,
    const StepCallback& stepCallback)
{
    // ---------- 功能块 0：输入与齐次初值检查 ----------
    const int dofs = numDofs(space);
    validateStateVector(initialDisplacements, dofs, "initial displacements");
    validateStateVector(initialVelocities, dofs, "initial velocities");
    requireHomogeneousConstraints(initialDisplacements, constrainedDofs,
                                  "initial displacements");
    requireHomogeneousConstraints(initialVelocities, constrainedDofs,
                                  "initial velocities");
    if (!std::isfinite(settings.timeStep) || settings.timeStep <= 0.0
        || settings.numberOfSteps < 0 || settings.maxIterations <= 0
        || settings.relativeTolerance < 0.0 || settings.absoluteTolerance < 0.0) {
        throw std::invalid_argument("invalid dynamics settings");
    }

    // ---------- 功能块 1（全程一次）：质量、参数与初始加速度 ----------
    DynamicsResult result;
    result.massMatrix = assembleConsistentMass(space, referenceDensity);
    result.massDiagnostics = diagnoseMassMatrix(space, referenceDensity,
                                                result.massMatrix);
    if (result.massDiagnostics.relativeMassError > 1.0e-11
        || result.massDiagnostics.relativeSymmetryError > 1.0e-13
        || !result.massDiagnostics.positiveDefinite) {
        throw std::runtime_error("consistent mass matrix self-check failed");
    }
    const GeneralizedAlphaParameters parameters =
        generalizedAlphaFromSpectralRadius(settings.spectralRadius);
    Eigen::VectorXd acceleration = computeInitialAcceleration(
        space, material, result.massMatrix, initialDisplacements,
        constrainedDofs, externalLoad, settings.elementFormulation);
    // 计算初始时间步的状态
    DynamicState current = makeState(0.0, initialDisplacements, initialVelocities,
                                     acceleration, result.massMatrix, space, material);
    // states[0] 是 t=0 的物理状态；回调可以在任何时间步安全请求停止。
    result.states.reserve(static_cast<std::size_t>(settings.numberOfSteps + 1));
    result.states.push_back(current);
    if (stepCallback && !stepCallback(result.states.back())) {
        return result;
    }

    // ---------- 功能块 2：外层时间步循环 ----------
    const double dt = settings.timeStep;
    const double dt2 = dt * dt;
    for (int step = 1; step <= settings.numberOfSteps; ++step) {
        // 同速度预测子：令 a_{n+1}^{(0)}=(gamma-1)/gamma*a_n，并直接取
        // v_{n+1}^{(0)}=v_n；d_{n+1}^{(0)} 仍由 Newmark 位移式确定。
        // &预测器predictor：同速度预测子
        Eigen::VectorXd nextAcceleration =
            ((parameters.gamma - 1.0) / parameters.gamma) * current.accelerations;
        Eigen::VectorXd nextDisplacement = current.displacements
            + dt * current.velocities
            + dt2 * ((0.5 - parameters.beta) * current.accelerations
                     + parameters.beta * nextAcceleration);
        Eigen::VectorXd nextVelocity = current.velocities;

        // ---------- 功能块 3：当前时间步的 Newton 校正 ----------
        // 未知量始终是 a_{n+1}。d_{n+1}、v_{n+1} 与它同步更新，
        // 但残量只通过 evaluateDynamicResidual() 在中间构形上求值。
        // &多步修正corrector
        bool converged = false;
        double initialResidualNorm = 0.0;
        for (int iteration = 1; iteration <= settings.maxIterations; ++iteration) {
            const DynamicResidualResult dynamic = evaluateDynamicResidual(
                space, material, result.massMatrix, current, nextAcceleration,
                dt, parameters, externalLoad, settings.elementFormulation);// &M^* R
            const double residualNorm = freeDofNorm(dynamic.residual, constrainedDofs); //计算自由度残差
            if (iteration == 1) {
                initialResidualNorm = residualNorm;
            }
            const double relativeResidual = initialResidualNorm > 0.0
                ? residualNorm / initialResidualNorm : 0.0;
            result.history.push_back({step, iteration, current.time + dt,
                                      residualNorm, relativeResidual});
            if (settings.printIterations) {
                std::cout << "step " << step << ", iteration " << iteration
                          << ": |R|=" << residualNorm
                          << ", relative=" << relativeResidual << '\n';
            }
            // 动力学可能没有外载荷，故不能用 ||f_ext|| 作分母；采用本时间步
            // 预测残量 R_0 的相对范数，与绝对范数取“或”。
            if (residualNorm < settings.absoluteTolerance
                || (iteration > 1
                    && relativeResidual < settings.relativeTolerance)) {
                converged = true;
                break;
            }
            // &施加dirichlet边界条件，开始牛顿迭代
            const ConstrainedLinearSystem constrained = applyHomogeneousDirichlet(
                dynamic.effectiveTangent, -dynamic.residual, constrainedDofs);
            // M* 对称但非线性切线未必全程正定，因此使用 SimplicialLDLT。
            Eigen::SimplicialLDLT<SparseMatrix> solver;
            solver.compute(constrained.matrix);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("effective tangent factorization failed at step "
                                         + std::to_string(step));
            }
            const Eigen::VectorXd correction = solver.solve(constrained.rightHandSide);
            if (solver.info() != Eigen::Success || !correction.allFinite()) {
                throw std::runtime_error("effective tangent solve failed at step "
                                         + std::to_string(step));
            }
            // 同一个 Delta a 同步校正三个 n+1 状态；这里没有 alpha_f。
            nextAcceleration += correction;
            nextDisplacement += parameters.beta * dt2 * correction;
            nextVelocity += parameters.gamma * dt * correction;
            // admissibility 检查使用本次残量对应的中间构形，而不是渲染状态。
            const Eigen::VectorXd displacementAlphaF = current.displacements
                + parameters.alphaF * (nextDisplacement - current.displacements);
            if (!hasPositiveDeformationJacobians(space, displacementAlphaF)) {
                throw std::runtime_error("non-positive deformation Jacobian at step "
                                         + std::to_string(step));
            }
        }
        if (!converged) {
            throw std::runtime_error("generalized-alpha Newton iteration did not converge at step "
                                     + std::to_string(step));
        }
        // ---------- 功能块 4：提交已收敛的物理状态 ----------
        // 只有走到这里，n+1 才成为新的 n；Newton 中间状态绝不进入历史。
        current = makeState(current.time + dt, std::move(nextDisplacement),
                            std::move(nextVelocity), std::move(nextAcceleration),
                            result.massMatrix, space, material);
        result.states.push_back(current);
        if (stepCallback && !stepCallback(result.states.back())) {
            break;
        }
    }
    return result;
}

double computeStrainEnergy(const nurbs::TensorProductNurbsSpace& space,
                           const Eigen::VectorXd& displacements,
                           const HyperelasticMaterial& material)
{
    validateStateVector(displacements, numDofs(space), "displacements");
    // ---------- 功能块 1：在参考域逐单元、逐 Gauss 点积分材料应变能 ----------
    // Psi 由具体材料直接提供；不能再用只对线性 SVK 成立的 0.5*S:E 代替。
    double energy = 0.0;
    for (int element = 0; element < space.numElements(); ++element) {
        const nurbs::ControlPointMatrix local = gatherElementDisplacements(
            space, element, displacements);
        for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
            for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
                for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                    const ReferencePointData reference = evaluateReferencePoint(
                        space, element, {q0, q1, q2});
                    const KinematicState state = evaluateKinematics(
                        local, reference.referenceGradients);
                    const MaterialResponse response = material.evaluate(state);
                    if (!std::isfinite(response.strainEnergyDensity)) {
                        throw std::domain_error(
                            "computeStrainEnergy: material energy density must be finite");
                    }
                    energy += reference.integrationWeight
                        * response.strainEnergyDensity;
                }
            }
        }
    }
    return energy;
}

} // namespace iga
