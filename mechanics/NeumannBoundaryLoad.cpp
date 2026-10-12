#include "mechanics/NeumannBoundaryLoad.h"

#include <Eigen/Geometry>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace iga {
namespace {

double mapGaussPointToSpan(double parentPoint, double lower, double upper)
{
    return 0.5 * ((upper - lower) * parentPoint + upper + lower);
}

double parameterMinimum(const nurbs::TensorProductNurbsSpace& space, int direction)
{
    const int firstSpan = space.span(direction, 0);
    return space.knotVector(direction)[firstSpan];
}

double parameterMaximum(const nurbs::TensorProductNurbsSpace& space, int direction)
{
    const int lastElement = space.numElements(direction) - 1;
    const int lastSpan = space.span(direction, lastElement);
    return space.knotVector(direction)[lastSpan + 1];
}

} // namespace

TimeDependentLoad makeReferenceFaceTractionLoad(
    const nurbs::TensorProductNurbsSpace& space,
    int normalDirection,
    int side,
    ReferenceTraction traction)
{
    if (normalDirection < 0 || normalDirection >= 3) {
        throw std::invalid_argument(
            "makeReferenceFaceTractionLoad: normal direction must be 0, 1, or 2");
    }
    if (side != 0 && side != 1) {
        throw std::invalid_argument(
            "makeReferenceFaceTractionLoad: side must be 0 or 1");
    }
    if (!traction) {
        throw std::invalid_argument(
            "makeReferenceFaceTractionLoad: traction callback must not be empty");
    }

    // ---------- 功能块 1：确定参数面及其两个自由积分方向 ----------
    // normalDirection 只表示哪个参数保持常量，不等同于物理空间的法向轴。
    const int free0 = (normalDirection + 1) % 3;
    const int free1 = (normalDirection + 2) % 3;
    std::array<double, 3> parameters{};
    for (int direction = 0; direction < 3; ++direction) {
        parameters[direction] = 0.5
            * (parameterMinimum(space, direction)
               + parameterMaximum(space, direction));
    }
    parameters[normalDirection] = side == 0
        ? parameterMinimum(space, normalDirection)
        : parameterMaximum(space, normalDirection);

    // ---------- 功能块 2：一次性积分每个全局基函数的参考面积系数 ----------
    // 面内两个方向逐非零 knot span 积分；固定方向位于 patch 端点，所以使用
    // 任意参数点求值路径。这里的额外基函数计算只发生在载荷创建阶段。
    Eigen::VectorXd nodalAreaWeights = Eigen::VectorXd::Zero(
        space.controlPoints().rows());
    for (int element1 = 0; element1 < space.numElements(free1); ++element1) {
        const int span1 = space.span(free1, element1);
        const double lower1 = space.knotVector(free1)[span1];
        const double upper1 = space.knotVector(free1)[span1 + 1];
        const double parentJacobian1 = 0.5 * (upper1 - lower1);
        for (int element0 = 0; element0 < space.numElements(free0); ++element0) {
            const int span0 = space.span(free0, element0);
            const double lower0 = space.knotVector(free0)[span0];
            const double upper0 = space.knotVector(free0)[span0 + 1];
            const double parentJacobian0 = 0.5 * (upper0 - lower0);

            for (int q1 = 0; q1 < space.numGaussPoints(free1); ++q1) {
                parameters[free1] = mapGaussPointToSpan(
                    space.gaussPoint(free1, q1), lower1, upper1);
                for (int q0 = 0; q0 < space.numGaussPoints(free0); ++q0) {
                    parameters[free0] = mapGaussPointToSpan(
                        space.gaussPoint(free0, q0), lower0, upper0);
                    const nurbs::PointShapeFunctionData point =
                        space.evaluateShapeFunctionsAt(parameters);

                    // 功能块 2.1：由参考控制点和参数导数组成两个物理切向量。
                    nurbs::ControlPointMatrix localPoints(point.connectivity.size(), 3);
                    for (std::size_t local = 0; local < point.connectivity.size(); ++local) {
                        localPoints.row(static_cast<Eigen::Index>(local)) =
                            space.controlPoints().row(point.connectivity[local]);
                    }
                    const Eigen::Vector3d tangent0 = localPoints.transpose()
                        * point.shape.parametricGradients.col(free0);
                    const Eigen::Vector3d tangent1 = localPoints.transpose()
                        * point.shape.parametricGradients.col(free1);
                    const double surfaceJacobian = tangent0.cross(tangent1).norm();
                    if (!std::isfinite(surfaceJacobian) || surfaceJacobian <= 0.0) {
                        throw std::domain_error(
                            "makeReferenceFaceTractionLoad: reference face Jacobian must be positive");
                    }

                    // 功能块 2.2：完整权重包含参考面积 Jacobian、两个 span 映射
                    // Jacobian 以及两个一维 Gauss 权重。
                    const double integrationWeight = surfaceJacobian
                        * parentJacobian0 * parentJacobian1
                        * space.gaussWeight(free0, q0)
                        * space.gaussWeight(free1, q1);
                    for (Eigen::Index local = 0; local < point.shape.values.size(); ++local) {
                        const int global = point.connectivity[static_cast<std::size_t>(local)];
                        nodalAreaWeights(global) += point.shape.values(local)
                            * integrationWeight;
                    }
                }
            }
        }
    }
    if (!nodalAreaWeights.allFinite() || nodalAreaWeights.sum() <= 0.0) {
        throw std::domain_error(
            "makeReferenceFaceTractionLoad: integrated reference face area is invalid");
    }

    // ---------- 功能块 3：返回只执行时间求值和自由度展开的轻量回调 ----------
    // 回调只捕获预积分系数，不引用 NURBS space，因此其生命周期与 space 无关。
    return [weights = std::move(nodalAreaWeights),
            traction = std::move(traction)](double time) {
        if (!std::isfinite(time)) {
            throw std::invalid_argument(
                "reference face traction load: time must be finite");
        }
        const Eigen::Vector3d tractionValue = traction(time);
        if (!tractionValue.allFinite()) {
            throw std::domain_error(
                "reference face traction load: traction callback returned non-finite values");
        }
        Eigen::VectorXd force = Eigen::VectorXd::Zero(3 * weights.size());
        for (Eigen::Index basis = 0; basis < weights.size(); ++basis) {
            force.segment<3>(3 * basis) = weights(basis) * tractionValue;
        }
        return force;
    };
}

} // namespace iga
