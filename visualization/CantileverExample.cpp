#include "visualization/CantileverExample.h"

#include "mechanics/GlobalNonlinearSystem.h"

#include "mechanics/SaintVenantKirchhoffMaterial.h"
#include "mechanics/NeoHookeanMaterial.h"
#include "mechanics/PrincipalStretchMaterial.h"

#include "visualization/BoundarySampler.h"

#include <array>
#include <utility>
#include <vector>

namespace visualization {
namespace {

// 生成 [0,1] 上的开放均匀 knot vector；端点重复 p+1 次。
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

// Greville 横坐标用于构造精确的线性长方体几何控制网格。
std::vector<double> grevilleAbscissae(const std::vector<double>& knots, int degree)
{
    const int count = static_cast<int>(knots.size()) - degree - 1;
    std::vector<double> values(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        for (int j = 1; j <= degree; ++j) {
            values[static_cast<std::size_t>(i)] += knots[static_cast<std::size_t>(i + j)];
        }
        values[static_cast<std::size_t>(i)] /= degree;
    }
    return values;
}

} // namespace

// 构造 10x1x1 的三次单 patch；xi 方向沿梁轴，zeta 反向为重力方向。
nurbs::TensorProductNurbsSpace makeCantileverSpace()
{
    const std::array<int, 3> degree{3, 3, 3};
    const std::array<int, 3> elements{16, 2, 2};
    const Eigen::Vector3d dimensions(10.0, 1.0, 1.0);
    std::array<std::vector<double>, 3> knots{
        openUniformKnots(degree[0], elements[0]),
        openUniformKnots(degree[1], elements[1]),
        openUniformKnots(degree[2], elements[2])};
    std::array<std::vector<double>, 3> greville{
        grevilleAbscissae(knots[0], degree[0]),
        grevilleAbscissae(knots[1], degree[1]),
        grevilleAbscissae(knots[2], degree[2])};
    const int n0 = static_cast<int>(greville[0].size());
    const int n1 = static_cast<int>(greville[1].size());
    const int n2 = static_cast<int>(greville[2].size());
    // 控制点全局编号与求解器一致：xi 最快，其次 eta，最后 zeta。
    nurbs::ControlPointMatrix points(n0 * n1 * n2, 3);
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
    return nurbs::TensorProductNurbsSpace(
        {nurbs::KnotVector(std::move(knots[0])),
         nurbs::KnotVector(std::move(knots[1])),
         nurbs::KnotVector(std::move(knots[2]))},
        degree, std::move(points), Eigen::VectorXd::Ones(global));
}

CantileverVisualization solveCantileverVisualization(int samplesPerSpan)
{
    // 功能块 1：建立几何、SVK 材料和沿 -z 的一致体力载荷。
    const nurbs::TensorProductNurbsSpace space = makeCantileverSpace();

    // const iga::SaintVenantKirchhoffMaterial material(0.0, 5.0e6);
    // const iga::NeoHookeanMaterial material(0.0, 5.0e6);
    const iga::PrincipalStretchMaterial material(0.0, 5.0e6);

    const Eigen::Vector3d gravity(0.0, 0.0, -1.0);

    // rho=500 is the already-validated nonlinear self-weight case from the global tests.
    const Eigen::VectorXd force = iga::assembleBodyForce(space, 500.0, gravity);

    // 功能块 2：固定 xi-min 端，通过现有 Newton 求解器得到控制点位移。
    iga::NewtonSettings settings;
    settings.printIterations = true;
    settings.maxIterations = 30;
    settings.elementFormulation = iga::ElementFormulation::UpdatedLagrangian;
    const iga::NewtonResult solution = iga::solveDeadLoad(
        space, material, force, iga::collectXiMinFaceDofs(space), settings);

    // 功能块 3：真实位移进入应力求值；轴向明确指定为全局 x 方向。
    FieldEvaluator evaluator(space, solution.displacements, material,
                             Eigen::Vector3d::UnitX(), gravity);
    CantileverVisualization result;
    // 后处理采样密度独立于上面的分析网格和求解精度。
    result.mesh = BoundarySampler(samplesPerSpan).sample(space, evaluator);
    // 在自由端截面中心记录位移，便于界面和测试快速检查重力方向是否正确。
    result.tipDisplacement = iga::evaluateDisplacement(
        space, solution.displacements, {1.0, 0.5, 0.5})(2);
    result.newtonIterations = static_cast<int>(solution.history.size());
    return result;
}

} // namespace visualization
