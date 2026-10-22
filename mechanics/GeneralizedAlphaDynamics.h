#ifndef GENERALIZEDALPHADYNAMICS_H
#define GENERALIZEDALPHADYNAMICS_H

#include "mechanics/ExternalLoad.h"
#include "mechanics/GlobalNonlinearSystem.h"

#include <functional>
#include <vector>

namespace iga {

struct GeneralizedAlphaParameters {
    /// 无限频率处的谱半径 rho_inf；越小，高频数值阻尼越强。
    double spectralRadius = 0.8;
    /// 加速度中间时刻参数 alpha_m。
    double alphaM = 2.0 / 3.0;
    /// 位移和外载荷中间时刻参数 alpha_f。
    double alphaF = 5.0 / 9.0;
    /// Newmark 速度更新参数，由 alpha_m、alpha_f 的二阶条件确定。
    double gamma = 11.0 / 18.0;
    /// Newmark 位移更新参数，由 alpha_m、alpha_f 的二阶条件确定。
    double beta = 25.0 / 81.0;
};

/**
 * @brief 由 rho_inf 生成满足二阶精度和无条件稳定条件的 generalized-alpha 参数。
 * @param spectralRadius 无限频率谱半径，范围为 [0,1]。
 * @return 同一套 alpha 约定下的 alpha_m、alpha_f、gamma 和 beta。
 */
GeneralizedAlphaParameters generalizedAlphaFromSpectralRadius(double spectralRadius);

/**
 * @brief 在参考构形中组装完整一致质量矩阵 M，不进行质量集中化。
 * @param space 三维张量积 NURBS 空间及参考几何。
 * @param referenceDensity 参考密度 rho_0，必须为正数。
 * @return N_dof x N_dof 稀疏矩阵；每个控制点对块为标量质量乘 I_3。
 *
 * M 只依赖参考构形和 rho_0，应在整个时间积分开始前计算一次。
 */
SparseMatrix assembleConsistentMass(const nurbs::TensorProductNurbsSpace& space,
                                    double referenceDensity);

/** @brief 一致质量矩阵的常驻自检结果。 */
struct MassDiagnostics {
    /// 由 1^T M 1 / 3 得到的离散总质量。
    double integratedMass = 0.0;
    /// 由参考域体积分独立计算的 rho_0 V_0。
    double expectedMass = 0.0;
    /// integratedMass 与 expectedMass 的相对误差。
    double relativeMassError = 0.0;
    /// ||M-M^T||/||M||，用于发现只填充半个对称矩阵的问题。
    double relativeSymmetryError = 0.0;
    /// SimplicialLLT 是否成功；完整一致质量矩阵应为对称正定。
    bool positiveDefinite = false;
};

/**
 * @brief 检查一致质量矩阵的总质量、对称性和正定性。
 * @return 诊断量；本函数不修改传入矩阵，也不施加边界条件。
 */
MassDiagnostics diagnoseMassMatrix(const nurbs::TensorProductNurbsSpace& space,
                                   double referenceDensity,
                                   const SparseMatrix& massMatrix);

/** @brief 一个已收敛物理时刻的完整状态，而非 n+alpha_f 中间构形。 */
struct DynamicState {
    /// 当前物理时间 t_n。
    double time = 0.0;
    /// d_n、v_n、a_n 均为 N_dof x 1，采用 I=3*A+i 的交错布局。
    Eigen::VectorXd displacements;
    Eigen::VectorXd velocities;
    Eigen::VectorXd accelerations;
    /// 0.5*v^T*M*v 与积分 Psi，供能量诊断和可视化使用。
    double kineticEnergy = 0.0;
    double strainEnergy = 0.0;
};

/** @brief 单个时间步内一次 Newton 迭代的收敛记录。 */
struct DynamicIterationRecord {
    int timeStep = 0;
    int iteration = 0;
    double time = 0.0;
    double residualNorm = 0.0;
    double relativeResidual = 0.0;
};

/** @brief generalized-alpha 时间循环和每步 Newton 迭代配置。 */
struct DynamicsSettings {
    /// 固定时间步长 Delta t。
    double timeStep = 1.0e-3;
    /// 从初始状态开始推进的时间步数，不包含第 0 帧。
    int numberOfSteps = 1;
    /// 每个时间步允许的最大 Newton 次数。
    int maxIterations = 20;
    /// ||R_i||/||R_0|| 和 ||R_i|| 的“或”收敛阈值。
    double relativeTolerance = 1.0e-8;
    double absoluteTolerance = 1.0e-12;
    /// rho_inf；默认 0.8 用于缓慢滤除阶跃载荷激发的高频模态。
    double spectralRadius = 0.8;
    bool printIterations = true;
    /// 单元内力和切线的积分构形；质量矩阵仍只在参考构形装配一次。
    ElementFormulation elementFormulation = ElementFormulation::TotalLagrangian;
};

/// 每次提交已收敛状态后调用；返回 false 可安全提前结束时间循环。
using StepCallback = std::function<bool(const DynamicState&)>;

/** @brief 动力学求解结果：常量质量矩阵、全部物理状态和 Newton 历史。 */
struct DynamicsResult {
    SparseMatrix massMatrix;
    MassDiagnostics massDiagnostics;
    std::vector<DynamicState> states;
    std::vector<DynamicIterationRecord> history;
};

/** @brief 固定 a_{n+1} 时的离散残量及其对 a_{n+1} 的精确导数。 */
struct DynamicResidualResult {
    /// R=M*a_{n+alpha_m}+f_int(d_{n+alpha_f})-f_ext(t_{n+alpha_f})。
    Eigen::VectorXd residual;
    /// M*=alpha_m*M+alpha_f*beta*Delta_t^2*K_T。
    SparseMatrix effectiveTangent;
    /// 本次残量实际使用的中间位移和中间加速度，主要用于测试和诊断。
    Eigen::VectorXd displacementAlphaF;
    Eigen::VectorXd accelerationAlphaM;
};

/**
 * @brief 把 a_{n+1} 视为 Newton 未知量，计算 R(a_{n+1}) 和有效切线 M*。
 *
 * 该函数复用静力层 assembleInternalSystem()，但求值位移是
 * d_{n+alpha_f}；外载荷也严格取在 t_{n+alpha_f}。
 */
DynamicResidualResult evaluateDynamicResidual(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const SparseMatrix& massMatrix,
    const DynamicState& previousState,
    const Eigen::VectorXd& nextAcceleration,
    double timeStep,
    const GeneralizedAlphaParameters& parameters,
    const TimeDependentLoad& externalLoad,
    ElementFormulation formulation = ElementFormulation::TotalLagrangian);

/**
 * @brief 在自由自由度上求解 M*a_0=f_ext(0)-f_int(d_0)。
 *
 * 受约束自由度通过齐次 Dirichlet 处理保持 a_0=0；除非初态已经静力
 * 平衡，不能直接把 a_0 设为零。
 */
Eigen::VectorXd computeInitialAcceleration(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const SparseMatrix& massMatrix,
    const Eigen::VectorXd& initialDisplacements,
    const std::vector<int>& constrainedDofs,
    const TimeDependentLoad& externalLoad,
    ElementFormulation formulation = ElementFormulation::TotalLagrangian);

/**
 * @brief 用 generalized-alpha predictor-multicorrector 推进非线性结构动力学。
 *
 * 函数内部只组装一次 M；每个时间步以 a_{n+1} 为 Newton 未知量，只有
 * 收敛后的 (d_{n+1},v_{n+1},a_{n+1}) 才写入 states 并传给回调。
 */
DynamicsResult solveDynamics(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    double referenceDensity,
    const Eigen::VectorXd& initialDisplacements,
    const Eigen::VectorXd& initialVelocities,
    const std::vector<int>& constrainedDofs,
    const TimeDependentLoad& externalLoad,
    const DynamicsSettings& settings = DynamicsSettings{},
    const StepCallback& stepCallback = StepCallback{});

/**
 * @brief 在参考域积分材料模型给出的超弹性应变能密度 Psi。
 * @return 当前已收敛位移状态的总应变能。
 */
double computeStrainEnergy(const nurbs::TensorProductNurbsSpace& space,
                           const Eigen::VectorXd& displacements,
                           const HyperelasticMaterial& material);

} // namespace iga

#endif // GENERALIZEDALPHADYNAMICS_H
