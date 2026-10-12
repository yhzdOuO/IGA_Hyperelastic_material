#ifndef GLOBALNONLINEARSYSTEM_H
#define GLOBALNONLINEARSYSTEM_H

#include "mechanics/TotalLagrangianElement.h"

#include <Eigen/SparseCore>
#include <array>
#include <vector>

namespace iga {

using SparseMatrix = Eigen::SparseMatrix<double>;

struct GlobalAssemblyResult {
    /// 当前全局位移下的 N_dof x N_dof 稀疏一致切线 K_T。
    SparseMatrix tangentStiffness;
    /// 当前全局位移下的 N_dof x 1 内力 f_int。
    Eigen::VectorXd internalForce;
};

struct ConstrainedLinearSystem {
    /// 约束行列清零并补单位对角元后的对称矩阵。
    SparseMatrix matrix;
    /// 约束位置置零后的右端项。
    Eigen::VectorXd rightHandSide;
};

/** @brief 静力载荷步进、Newton 迭代和可接受性线搜索配置。 */
struct NewtonSettings {
    int loadSteps = 1;
    int maxIterations = 25;
    int maxLineSearchCuts = 12;
    double relativeTolerance = 1.0e-9;
    double absoluteTolerance = 1.0e-12;
    bool printIterations = true;
};

/** @brief 静力求解中一次 Newton 迭代的残量记录。 */
struct NewtonIterationRecord {
    int loadStep = 0;
    int iteration = 0;
    double residualNorm = 0.0;
    double relativeResidual = 0.0;
};

/** @brief 静力求解最终位移、未约束残量及完整迭代历史。 */
struct NewtonResult {
    Eigen::VectorXd displacements;
    Eigen::VectorXd residual;
    std::vector<NewtonIterationRecord> history;
};

/**
 * @brief 返回全局位移自由度数 N_dof。
 * @param space 待计数的 NURBS patch。
 * @return 3*N_cp；自由度布局固定为 I=3*A+i。
 */
int numDofs(const nurbs::TensorProductNurbsSpace& space);

/**
 * @brief 按 I=3*A+i 从全局向量收集一个单元的控制变量。
 * @param space 提供惰性 IEN 连接关系的 NURBS patch。
 * @param element 单元编号。
 * @param globalDisplacements N_dof x 1 全局交错位移向量。
 * @return n_en x 3 局部位移矩阵，每行对应一个活动基函数。
 */
nurbs::ControlPointMatrix gatherElementDisplacements(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const Eigen::VectorXd& globalDisplacements);

/**
 * @brief 在参考构形中组装一致体力向量。
 * @param space NURBS patch 及参考几何。
 * @param referenceDensity 参考密度 rho_0。
 * @param acceleration 常量死载荷加速度 g。
 * @return integral(R_A*rho_0*g)dOmega_0 组成的 N_dof x 1 向量。
 */
Eigen::VectorXd assembleBodyForce(
    const nurbs::TensorProductNurbsSpace& space,
    double referenceDensity,
    const Eigen::Vector3d& acceleration);

/**
 * @brief 用 triplet 装配当前构形的全局切线和内力。
 * @param space 提供几何和连接关系的 NURBS patch。
 * @param globalDisplacements 当前交错位移向量，I=3*A+i。
 * @param material 所有单元共用的材料模型。
 * @return 未施加边界条件的稀疏 K_T 与 f_int；静力和动力学共同复用。
 */
GlobalAssemblyResult assembleInternalSystem(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements,
    const HyperelasticMaterial& material);

/**
 * @brief 收集 xi 最小控制层上三个方向的全部位移自由度。
 * @param space 使用开放结点向量的 NURBS patch。
 * @return 按 I=3*A+i 排列的约束自由度列表。
 */
std::vector<int> collectXiMinFaceDofs(const nurbs::TensorProductNurbsSpace& space);

/**
 * @brief 在不破坏对称性的前提下施加齐次 Dirichlet 增量条件。
 * @param matrix 未约束的稀疏矩阵，可为 K_T、M 或有效切线 M*。
 * @param rightHandSide 未约束线性系统右端项。
 * @param constrainedDofs 增量必须保持为零的自由度。
 * @return 约束行列清零、对角置 1 且相应 RHS 置零的线性系统。
 */
ConstrainedLinearSystem applyHomogeneousDirichlet(
    const SparseMatrix& matrix,
    const Eigen::VectorXd& rightHandSide,
    const std::vector<int>& constrainedDofs);

/**
 * @brief 只在自由自由度上计算全局向量的欧氏范数。
 * @param vector N_dof x 1 全局向量。
 * @param constrainedDofs 从范数中排除的自由度；其残量通常是支反力。
 */
double freeDofNorm(const Eigen::VectorXd& vector,
                   const std::vector<int>& constrainedDofs);

/**
 * @brief 检查试探位移下每个单元 Gauss 点是否满足 det(F)>0。
 * @return 仅当所有变形梯度有限且保持局部定向时返回 true。
 */
bool hasPositiveDeformationJacobians(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements);

/**
 * @brief 在任意参数点计算物理位移场，而不是读取某个控制变量。
 * @param parameters 参数坐标 {xi,eta,zeta}。
 * @return sum_A R_A*d_A 得到的三维插值位移。
 */
Eigen::Vector3d evaluateDisplacement(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements,
    const std::array<double, 3>& parameters);

/**
 * @brief 用载荷分级 Newton 法求解齐次 Dirichlet 静力死载荷问题。
 * @param totalExternalForce 迭代前一次组装的完整死载荷向量。
 * @param settings 载荷分级、收敛阈值和可接受性线搜索设置。
 * @return 收敛位移、最终未约束残量和迭代历史。
 */
NewtonResult solveDeadLoad(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const Eigen::VectorXd& totalExternalForce,
    const std::vector<int>& constrainedDofs,
    const NewtonSettings& settings = NewtonSettings{});

} // namespace iga

#endif // GLOBALNONLINEARSYSTEM_H
