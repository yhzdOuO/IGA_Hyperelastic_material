#include "mechanics/GlobalNonlinearSystem.h"

#include <Eigen/SparseCholesky>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace iga {

namespace {

std::vector<unsigned char> constrainedMask(
    int count, const std::vector<int>& constrainedDofs)
{
    // 用长度为 N_dof 的字节数组标记约束，后续判断为 O(1)。
    std::vector<unsigned char> mask(static_cast<std::size_t>(count), 0);
    for (const int dof : constrainedDofs) {
        if (dof < 0 || dof >= count) {
            throw std::out_of_range("constrained DOF is out of range");
        }
        mask[static_cast<std::size_t>(dof)] = 1;
    }
    return mask; //&类似[0, 0, 1, 0, 1, 0]
}

void validateGlobalVector(const nurbs::TensorProductNurbsSpace& space,
                          const Eigen::VectorXd& vector,
                          const char* message)
{
    if (vector.size() != numDofs(space) || !vector.allFinite()) {
        throw std::invalid_argument(message);
    }
}

} // namespace

int numDofs(const nurbs::TensorProductNurbsSpace& space)
{
    // 每个控制变量含 x、y、z 三个分量，全局编号固定为 I=3*A+i。
    return 3 * static_cast<int>(space.controlPoints().rows());
}

nurbs::ControlPointMatrix gatherElementDisplacements(
    const nurbs::TensorProductNurbsSpace& space,
    int element,
    const Eigen::VectorXd& globalDisplacements)
{
    validateGlobalVector(space, globalDisplacements,
                         "gatherElementDisplacements: global vector has invalid size or values");
    // 从全局向量 d(N_dof x 1)中收集出 u_e(n_en x 3)。
    nurbs::ControlPointMatrix local(space.numElementNodes(), 3);
    const std::vector<int> connectivity = space.elementConnectivity(element);
    for (int a = 0; a < space.numElementNodes(); ++a) {
        const int globalBasis = connectivity[static_cast<std::size_t>(a)];
        local.row(a) = globalDisplacements.segment<3>(3 * globalBasis).transpose();
    }
    return local;
}

Eigen::VectorXd assembleBodyForce(
    const nurbs::TensorProductNurbsSpace& space,
    double referenceDensity,
    const Eigen::Vector3d& acceleration)
{
    if (!std::isfinite(referenceDensity) || referenceDensity < 0.0
        || !acceleration.allFinite()) {
        throw std::invalid_argument("assembleBodyForce: density and acceleration must be finite");
    }
    // ---------- 一致体力向量：f_ext,A = integral(R_A*rho0*g)dV ----------
    // 重力是 TL 参考构形中的死载荷，结果为 N_dof x 1，只需在 Newton 之前计算一次。
    Eigen::VectorXd externalForce = Eigen::VectorXd::Zero(numDofs(space));
    const Eigen::Vector3d bodyForce = referenceDensity * acceleration;
    for (int element = 0; element < space.numElements(); ++element) {
        const std::vector<int> connectivity = space.elementConnectivity(element);
        for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
            for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
                for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                    const ReferencePointData point =
                        evaluateReferencePoint(space, element, {q0, q1, q2});
                    // 每个 R_A 是标量，乘 3x1 的 rho0*g 后散射到全局的三个连续自由度。
                    for (int a = 0; a < space.numElementNodes(); ++a) {
                        const int globalBasis = connectivity[static_cast<std::size_t>(a)];
                        externalForce.segment<3>(3 * globalBasis).noalias() +=
                            point.integrationWeight * point.shapeValues(a) * bodyForce;
                    }
                }
            }
        }
    }
    return externalForce;
}

GlobalAssemblyResult assembleInternalSystem(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements,
    const HyperelasticMaterial& material)
{
    validateGlobalVector(space, globalDisplacements,
                         "assembleInternalSystem: global displacement vector is invalid");
    // ---------- 当前位移下的全局切线与内力 ----------
    const int dofs = numDofs(space);
    const int localDofs = 3 * space.numElementNodes();
    std::vector<Eigen::Triplet<double>> triplets;
    // 每个 Ke 是 localDofs x localDofs。验证阶段保留全部条目，
    // 让 setFromTriplets 自动累加相邻单元贡献，不预计算 CSR 索引。
    triplets.reserve(static_cast<std::size_t>(space.numElements())
                     * static_cast<std::size_t>(localDofs)
                     * static_cast<std::size_t>(localDofs));

    // K 为 N_dof x N_dof 稀疏矩阵，f_int 为 N_dof x 1。
    GlobalAssemblyResult result;
    result.internalForce = Eigen::VectorXd::Zero(dofs);
    for (int element = 0; element < space.numElements(); ++element) {
        const nurbs::ControlPointMatrix localDisplacements =
            gatherElementDisplacements(space, element, globalDisplacements);
        const ElementResult local =
            computeElementResponse(space, element, localDisplacements, material); // &计算一个element的切线刚度矩阵和内力向量
        const std::vector<int> connectivity = space.elementConnectivity(element); // &计算映射connectivity[local_index] = global_index，其实是IEN(e, a) = A省略了e
        // 用同一套 IEN 同时散射 Ke 的行列和 fe，避免自由度布局不一致。
        for (int a = 0; a < space.numElementNodes(); ++a) {
            const int globalBasisA = connectivity[static_cast<std::size_t>(a)];
            for (int componentA = 0; componentA < 3; ++componentA) {
                // 局部和全局均采用分量交错布局：local=3*a+i，global=3*A+i。
                const int globalRow = 3 * globalBasisA + componentA;
                const int localRow = 3 * a + componentA;
                result.internalForce(globalRow) += local.internalForce(localRow);
                for (int b = 0; b < space.numElementNodes(); ++b) {
                    const int globalBasisB = connectivity[static_cast<std::size_t>(b)];
                    for (int componentB = 0; componentB < 3; ++componentB) {
                        const int globalColumn = 3 * globalBasisB + componentB;
                        const int localColumn = 3 * b + componentB;
                        triplets.emplace_back(globalRow, globalColumn,
                                              local.tangentStiffness(localRow, localColumn));
                    }
                }
            }
        }
    }
    // 重复的 (I,J) triplet 会在这里自动求和，这正是有限元装配语义。
    result.tangentStiffness.resize(dofs, dofs);
    result.tangentStiffness.setFromTriplets(triplets.begin(), triplets.end());
    result.tangentStiffness.makeCompressed();
    return result; // &返回全局的切线刚度矩阵+内力向量
}

std::vector<int> collectXiMinFaceDofs(const nurbs::TensorProductNurbsSpace& space)
{
    const int n0 = space.numControlPoints(0);
    const int n1 = space.numControlPoints(1);
    const int n2 = space.numControlPoints(2);
    // xi=0 对应 i_xi=0 的整层控制点；开结点向量下可精确施加齐次位移。
    std::vector<int> result;
    result.reserve(static_cast<std::size_t>(3 * n1 * n2));
    for (int i2 = 0; i2 < n2; ++i2) {
        for (int i1 = 0; i1 < n1; ++i1) {
            const int globalBasis = i1 * n0 + i2 * n0 * n1;
            for (int component = 0; component < 3; ++component) {
                result.push_back(3 * globalBasis + component);
            }
        }
    }
    return result;
}

ConstrainedLinearSystem applyHomogeneousDirichlet(
    const SparseMatrix& matrix,
    const Eigen::VectorXd& rightHandSide,
    const std::vector<int>& constrainedDofs)
{
    if (matrix.rows() != matrix.cols() || matrix.rows() != rightHandSide.size()) {
        throw std::invalid_argument("applyHomogeneousDirichlet: system dimensions are inconsistent");
    }
    // ---------- 对 Newton 增量施加齐次 Dirichlet 条件 ----------
    // 跳过受约束 DOF 涉及的所有行列条目，再为它们补单位对角元。
    // 这样得到的矩阵仍为 N_dof x N_dof，并且保持对称。
    const int dofs = static_cast<int>(matrix.rows());
    const std::vector<unsigned char> mask = constrainedMask(dofs, constrainedDofs); //&类似[0, 0, 1, 0, 1, 0]
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(matrix.nonZeros()) + constrainedDofs.size());
    for (int outer = 0; outer < matrix.outerSize(); ++outer) {
        for (SparseMatrix::InnerIterator entry(matrix, outer); entry; ++entry) {
            if (mask[static_cast<std::size_t>(entry.row())] == 0
                && mask[static_cast<std::size_t>(entry.col())] == 0) {
                triplets.emplace_back(entry.row(), entry.col(), entry.value());
            } // &(i,j)中i或者j为限制自由度的节点时，M(i,j)=0，M(i,i)=1
        }
    }
    for (int dof = 0; dof < dofs; ++dof) {
        if (mask[static_cast<std::size_t>(dof)] != 0) {
            triplets.emplace_back(dof, dof, 1.0); // &M(i,i)=1
        }
    }

    ConstrainedLinearSystem result;
    result.matrix.resize(dofs, dofs);
    result.matrix.setFromTriplets(triplets.begin(), triplets.end());
    result.matrix.makeCompressed();
    // 受约束处令 rhs=0，保证求得的位移增量 delta_d(I)=0。
    result.rightHandSide = rightHandSide;
    for (int dof = 0; dof < dofs; ++dof) {
        if (mask[static_cast<std::size_t>(dof)] != 0) {
            result.rightHandSide(dof) = 0.0; // &R(i)=0
        }
    }
    return result;
}

double freeDofNorm(const Eigen::VectorXd& vector,
                   const std::vector<int>& constrainedDofs)
{
    if (!vector.allFinite()) {
        throw std::invalid_argument("freeDofNorm: vector must be finite");
    }
    const std::vector<unsigned char> mask =
        constrainedMask(static_cast<int>(vector.size()), constrainedDofs);
    // 约束自由度上的残量是支反力，不能计入 Newton 收敛范数。
    double squaredNorm = 0.0;
    for (int dof = 0; dof < vector.size(); ++dof) {
        if (mask[static_cast<std::size_t>(dof)] == 0) {
            squaredNorm += vector(dof) * vector(dof);
        }
    }
    return std::sqrt(squaredNorm);
}

bool hasPositiveDeformationJacobians(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements)
{
    validateGlobalVector(space, globalDisplacements,
                         "hasPositiveDeformationJacobians: displacement vector is invalid");
    // 在所有单元的所有 Gauss 点检查 F(3x3)=I+u_e^T*G，防止试探步翻转材料。
    for (int element = 0; element < space.numElements(); ++element) {
        const nurbs::ControlPointMatrix localDisplacements =
            gatherElementDisplacements(space, element, globalDisplacements);
        for (int q2 = 0; q2 < space.numGaussPoints(2); ++q2) {
            for (int q1 = 0; q1 < space.numGaussPoints(1); ++q1) {
                for (int q0 = 0; q0 < space.numGaussPoints(0); ++q0) {
                    const ReferencePointData point =
                        evaluateReferencePoint(space, element, {q0, q1, q2});
                    const Eigen::Matrix3d displacementGradient =
                        localDisplacements.transpose() * point.physicalGradients;
                    const Eigen::Matrix3d deformationGradient =
                        Eigen::Matrix3d::Identity() + displacementGradient;
                    const double determinant = deformationGradient.determinant();
                    if (!std::isfinite(determinant) || determinant <= 0.0) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

Eigen::Vector3d evaluateDisplacement(
    const nurbs::TensorProductNurbsSpace& space,
    const Eigen::VectorXd& globalDisplacements,
    const std::array<double, 3>& parameters)
{
    validateGlobalVector(space, globalDisplacements,
                         "evaluateDisplacement: global displacement vector is invalid");
    // 后处理必须计算 u(xi)=sum_A R_A(xi)d_A；控制变量本身不是物理节点值。
    const nurbs::PointShapeFunctionData point =
        space.evaluateShapeFunctionsAt(parameters);
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    for (int a = 0; a < point.shape.values.size(); ++a) {
        const int globalBasis = point.connectivity[static_cast<std::size_t>(a)];
        result.noalias() += point.shape.values(a)
                            * globalDisplacements.segment<3>(3 * globalBasis);
    }
    return result;
}

NewtonResult solveDeadLoad(
    const nurbs::TensorProductNurbsSpace& space,
    const HyperelasticMaterial& material,
    const Eigen::VectorXd& totalExternalForce,
    const std::vector<int>& constrainedDofs,
    const NewtonSettings& settings)
{
    validateGlobalVector(space, totalExternalForce,
                         "solveDeadLoad: external-force vector is invalid");
    if (settings.loadSteps <= 0 || settings.maxIterations <= 0
        || settings.maxLineSearchCuts < 0 || settings.relativeTolerance <= 0.0
        || settings.absoluteTolerance <= 0.0) {
        throw std::invalid_argument("solveDeadLoad: Newton settings are invalid");
    }
    constrainedMask(numDofs(space), constrainedDofs);

    // ---------- 初值与载荷分级 ----------
    // d 为 N_dof x 1，零初值天然满足本函数支持的齐次边界条件。
    NewtonResult result;
    result.displacements = Eigen::VectorXd::Zero(numDofs(space));
    result.residual = Eigen::VectorXd::Zero(numDofs(space));
    // 载荷步进：牛顿迭代之外的循环。载荷因子：逐步施加载荷，防止一开始载荷过大导致不收敛。 // & 载荷步+牛顿步（非线性静力学）—— 时间步+牛顿步（动力学）
    for (int loadStep = 1; loadStep <= settings.loadSteps; ++loadStep) {
        const double loadFactor = static_cast<double>(loadStep) / settings.loadSteps;
        const Eigen::VectorXd currentExternalForce = loadFactor * totalExternalForce;
        const double externalNorm = freeDofNorm(currentExternalForce, constrainedDofs);
        bool converged = false;
        // 每个载荷步以上一步收敛解为初值，并在 当前载荷因子 下做 Newton 迭代。
        for (int iteration = 1; iteration <= settings.maxIterations; ++iteration) {
            const GlobalAssemblyResult assembly = assembleInternalSystem(
                space, result.displacements, material);
            // 残量 R=f_int-f_ext，大小为 N_dof x 1；切线 K=dR/dd=d(f_int)/dd。
            result.residual = assembly.internalForce - currentExternalForce;
            const double residualNorm = freeDofNorm(result.residual, constrainedDofs);
            const double convergenceMeasure = externalNorm > 0.0
                ? residualNorm / externalNorm : residualNorm;
            result.history.push_back(
                {loadStep, iteration, residualNorm, convergenceMeasure});
            if (settings.printIterations) {
                std::cout << "load step " << loadStep << '/' << settings.loadSteps
                          << ", iter " << iteration
                          << ", rel_res = " << convergenceMeasure << '\n';
            }
            // 有外力时使用相对残量；零载荷问题退回绝对残量判据。
            const double tolerance = externalNorm > 0.0
                ? settings.relativeTolerance : settings.absoluteTolerance;
            if (convergenceMeasure < tolerance) {
                converged = true;
                break;
            }

            // 求解 K*delta_d=-R。SimplicialLDLT 可处理可能不定的对称切线，
            // 且每次 factorization/solve 后都必须检查 Eigen 返回状态。
            const ConstrainedLinearSystem constrained = applyHomogeneousDirichlet(
                assembly.tangentStiffness, -result.residual, constrainedDofs);
            Eigen::SimplicialLDLT<SparseMatrix> solver;
            solver.compute(constrained.matrix);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("solveDeadLoad: SimplicialLDLT factorization failed");
            }
            const Eigen::VectorXd increment = solver.solve(constrained.rightHandSide);
            if (solver.info() != Eigen::Success || !increment.allFinite()) {
                throw std::runtime_error("solveDeadLoad: SimplicialLDLT solve failed");
            }

            // ---------- 简易回退线搜索 ----------
            // 若完整 Newton 步导致任一积分点 det(F)<=0，则不断把步长减半。
            // F是当前构型x对参考构型X的雅可比，描述了变换的局部，它的行列式＞=0避免了局部翻转。
            double scale = 1.0;
            bool admissible = false;
            Eigen::VectorXd trial = result.displacements;
            for (int cut = 0; cut <= settings.maxLineSearchCuts; ++cut) {
                trial = result.displacements + scale * increment;
                if (hasPositiveDeformationJacobians(space, trial)) {
                    admissible = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!admissible) {
                throw std::runtime_error(
                    "solveDeadLoad: line search could not find a trial state with det(F)>0");
            }
            result.displacements = trial;
        }
        if (!converged) {
            throw std::runtime_error("solveDeadLoad: Newton iteration did not converge");
        }
    }
    return result;
}

} // namespace iga
