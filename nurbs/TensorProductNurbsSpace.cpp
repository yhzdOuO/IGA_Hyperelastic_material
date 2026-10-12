#include "nurbs/TensorProductNurbsSpace.h"

#include "nurbs/NurbsAlgorithm.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace nurbs {

namespace {

void requireDirection(int direction)
{
    if (direction < 0 || direction >= 3) {
        throw std::out_of_range("TensorProductNurbsSpace: direction must be 0, 1, or 2");
    }
}

void requireIndex(int index, int count, const char* message)
{
    if (index < 0 || index >= count) {
        throw std::out_of_range(message);
    }
}

} // namespace

TensorProductNurbsSpace::TensorProductNurbsSpace(
    std::array<KnotVector, 3> knotVectors,
    std::array<int, 3> degrees,
    ControlPointMatrix controlPoints,
    Eigen::VectorXd weights)
    : knotVectors_(std::move(knotVectors)),
      degrees_(degrees),
      controlPoints_(std::move(controlPoints)),
      weights_(std::move(weights))
{
    // 每个方向先由 m=n+p+1 推出控制点数量，再一次性建立 L0 数据。
    // 构造完成后，单元积分过程中不再搜索跨度或重复计算一维基函数。
    for (int direction = 0; direction < 3; ++direction) {
        controlPointCounts_[direction] =
            knotVectors_[direction].inferLastCtrlPtIndex(degrees_[direction]) + 1;
    }
    validate();
    buildSpans();
    buildTables();
}

const KnotVector& TensorProductNurbsSpace::knotVector(int direction) const
{
    requireDirection(direction);
    return knotVectors_[direction];
}

int TensorProductNurbsSpace::degree(int direction) const
{
    requireDirection(direction);
    return degrees_[direction];
}

int TensorProductNurbsSpace::numControlPoints(int direction) const
{
    requireDirection(direction);
    return controlPointCounts_[direction];
}

int TensorProductNurbsSpace::numElements(int direction) const
{
    requireDirection(direction);
    return static_cast<int>(spans_[direction].size());
}

int TensorProductNurbsSpace::numElements() const
{
    return numElements(0) * numElements(1) * numElements(2);
}

int TensorProductNurbsSpace::numElementNodes() const
{
    return (degrees_[0] + 1) * (degrees_[1] + 1) * (degrees_[2] + 1);
}

int TensorProductNurbsSpace::numGaussPoints(int direction) const
{
    requireDirection(direction);
    return static_cast<int>(gaussRules_[direction].points.size());
}

double TensorProductNurbsSpace::gaussPoint(int direction, int gauss) const
{
    requireDirection(direction);
    requireIndex(gauss, numGaussPoints(direction),
                 "TensorProductNurbsSpace: Gauss index is out of range");
    return gaussRules_[direction].points[static_cast<std::size_t>(gauss)];
}

double TensorProductNurbsSpace::gaussWeight(int direction, int gauss) const
{
    requireDirection(direction);
    requireIndex(gauss, numGaussPoints(direction),
                 "TensorProductNurbsSpace: Gauss index is out of range");
    return gaussRules_[direction].weights[static_cast<std::size_t>(gauss)];
}

int TensorProductNurbsSpace::span(int direction, int directionalElement) const
{
    requireDirection(direction);
    requireIndex(directionalElement, numElements(direction),
                 "TensorProductNurbsSpace: directional element is out of range");
    return spans_[direction][directionalElement];
}

std::array<int, 3> TensorProductNurbsSpace::elementIndices(int element) const
{
    requireIndex(element, numElements(), "TensorProductNurbsSpace: element is out of range");
    const int ne0 = numElements(0);
    const int ne1 = numElements(1);
    // 单元编号约定：xi 最快、eta 次之、zeta 最慢。
    std::array<int, 3> result{};
    result[0] = element % ne0;
    result[1] = (element / ne0) % ne1;
    result[2] = element / (ne0 * ne1);
    return result;
}

std::array<int, 3> TensorProductNurbsSpace::elementSpans(int element) const
{
    const std::array<int, 3> indices = elementIndices(element);
    return {spans_[0][indices[0]], spans_[1][indices[1]], spans_[2][indices[2]]};
}

int TensorProductNurbsSpace::IEN(int localBasis, int element) const
{
    requireIndex(localBasis, numElementNodes(),
                 "TensorProductNurbsSpace: local basis index is out of range");
    const std::array<int, 3> elementSpan = elementSpans(element);
    const int d0 = degrees_[0] + 1;
    const int d1 = degrees_[1] + 1;
    // 将局部扁平编号还原为 (a_xi,a_eta,a_zeta)，同样令 xi 最快变化。
    const int a0 = localBasis % d0;
    const int a1 = (localBasis / d0) % d1;
    const int a2 = localBasis / (d0 * d1);
    // base 是该单元局部 (0,0,0) 基函数对应的全局编号。
    const int base = (elementSpan[0] - degrees_[0])
                   + (elementSpan[1] - degrees_[1]) * controlPointCounts_[0]
                   + (elementSpan[2] - degrees_[2])
                     * controlPointCounts_[0] * controlPointCounts_[1];
    return base + a0 + a1 * controlPointCounts_[0]
           + a2 * controlPointCounts_[0] * controlPointCounts_[1];
}

std::vector<int> TensorProductNurbsSpace::elementConnectivity(int element) const
{
    const std::array<int, 3> elementSpan = elementSpans(element);
    const int base = (elementSpan[0] - degrees_[0])
                   + (elementSpan[1] - degrees_[1]) * controlPointCounts_[0]
                   + (elementSpan[2] - degrees_[2])
                     * controlPointCounts_[0] * controlPointCounts_[1];
    // 按 xi -> eta -> zeta 的顺序生成连接关系。内层只做 global++，
    // 避免对每一个局部基函数重复进行乘法和扁平化计算。
    std::vector<int> result(static_cast<std::size_t>(numElementNodes()));
    int local = 0;
    int planeStart = base;
    for (int a2 = 0; a2 <= degrees_[2]; ++a2) {
        int rowStart = planeStart;
        for (int a1 = 0; a1 <= degrees_[1]; ++a1) {
            int global = rowStart;
            for (int a0 = 0; a0 <= degrees_[0]; ++a0) {
                result[static_cast<std::size_t>(local++)] = global++;
            }
            rowStart += controlPointCounts_[0];
        }
        planeStart += controlPointCounts_[0] * controlPointCounts_[1];
    }
    return result;
}

ControlPointMatrix TensorProductNurbsSpace::elementControlPoints(int element) const
{
    const std::vector<int> connectivity = elementConnectivity(element);
    // result 大小为 n_en x 3，每行是一个局部控制点的参考坐标 X_A。
    ControlPointMatrix result(numElementNodes(), 3);
    for (int a = 0; a < numElementNodes(); ++a) {
        result.row(a) = controlPoints_.row(connectivity[static_cast<std::size_t>(a)]);
    }
    return result;
}

Eigen::VectorXd TensorProductNurbsSpace::elementWeights(int element) const
{
    const std::vector<int> connectivity = elementConnectivity(element);
    Eigen::VectorXd result(numElementNodes());
    for (int a = 0; a < numElementNodes(); ++a) {
        result(a) = weights_(connectivity[static_cast<std::size_t>(a)]);
    }
    return result;
}

ShapeFunctionData TensorProductNurbsSpace::evaluateShapeFunctions(
    int element, const std::array<int, 3>& gaussIndices) const
{
    const std::array<int, 3> elementIndex = elementIndices(element); // &从全局扁平的element index中找到三个方向上各自的element index
    for (int direction = 0; direction < 3; ++direction) {
        requireIndex(gaussIndices[direction], numGaussPoints(direction),
                     "TensorProductNurbsSpace: Gauss index is out of range");
    }

    // 本函数只查预计算表，不调用 FindSpan/DersBasisFuns。
    const Eigen::VectorXd localWeights = elementWeights(element);
    const int d0 = degrees_[0] + 1;
    const int d1 = degrees_[1] + 1;
    const int d2 = degrees_[2] + 1;
    // 先共享 eta-zeta 的乘积因子。三个数组长度均为 d_eta*d_zeta，
    // 可减少三维张量积中重复的乘法。
    const int bcCount = d1 * d2;
    std::vector<double> productBC(static_cast<std::size_t>(bcCount));
    std::vector<double> derivativeBC1(static_cast<std::size_t>(bcCount));
    std::vector<double> derivativeBC2(static_cast<std::size_t>(bcCount));

    for (int a2 = 0; a2 < d2; ++a2) {
        // &第2个方向的表中，找第2个方向的element index下，参数中高斯点的第2个方向的值对应的，第a2个基函数的0次导数
        const double value2 = tables_[2].value(elementIndex[2], gaussIndices[2], 0, a2);
        const double derivative2 = tables_[2].value(elementIndex[2], gaussIndices[2], 1, a2);
        for (int a1 = 0; a1 < d1; ++a1) {
            const int bc = a1 + a2 * d1;
            const double value1 = tables_[1].value(elementIndex[1], gaussIndices[1], 0, a1);
            const double derivative1 = tables_[1].value(elementIndex[1], gaussIndices[1], 1, a1);
            productBC[static_cast<std::size_t>(bc)] = value1 * value2;
            derivativeBC1[static_cast<std::size_t>(bc)] = derivative1 * value2; // &对M偏导 M'_{j}N_{k}
            derivativeBC2[static_cast<std::size_t>(bc)] = value1 * derivative2; // &对N偏导 M_{j}N'_{k}
        }
    }

    // weighted 为 n_en x 1 的加权 B 样条 W_A；weightedDerivatives 为
    // n_en x 3，每行依次存 dW_A/dxi、dW_A/deta、dW_A/dzeta。
    ShapeFunctionData result;
    result.values.resize(numElementNodes());
    result.parametricGradients.resize(numElementNodes(), 3);
    Eigen::VectorXd weighted(numElementNodes());
    Eigen::Matrix<double, Eigen::Dynamic, 3> weightedDerivatives(numElementNodes(), 3);
    Eigen::Vector3d derivativeSums = Eigen::Vector3d::Zero();
    double weightSum = 0.0;

    int local = 0;
    for (int a2 = 0; a2 < d2; ++a2) {
        for (int a1 = 0; a1 < d1; ++a1) {
            const int bc = a1 + a2 * d1;
            for (int a0 = 0; a0 < d0; ++a0) {
                const double value0 = tables_[0].value(elementIndex[0], gaussIndices[0], 0, a0);      // &L_{i}
                const double derivative0 = tables_[0].value(elementIndex[0], gaussIndices[0], 1, a0); // &L'_{i}
                const double weightedValue0 = localWeights(local) * value0; // & w_{ijk} L_{i}
                weighted(local) = weightedValue0 * productBC[static_cast<std::size_t>(bc)];    //& A_{l} = A_{ijk} = w_{ijk} L_{i} M_{j} N_{k}
                weightedDerivatives(local, 0) = localWeights(local) * derivative0
                                                * productBC[static_cast<std::size_t>(bc)];     //& W_{l} = W_{ijk} = w_{ijk} L'_{i} M_{j} N_{k}
                weightedDerivatives(local, 1) = weightedValue0
                                                * derivativeBC1[static_cast<std::size_t>(bc)]; //& B_{l} = B_{ijk} = w_{ijk} L_{i} M'_{j} N_{k}
                weightedDerivatives(local, 2) = weightedValue0
                                                * derivativeBC2[static_cast<std::size_t>(bc)]; //& C_{l} = C_{ijk} = w_{ijk} L_{i} M_{j} N'_{k}
                weightSum += weighted(local); // & S_{W}
                derivativeSums += weightedDerivatives.row(local).transpose(); // & [S_{A}, S_{B}, S_{C}]
                ++local;
            }
        }
    }

    // 有理化及商法则：R_A=W_A/Wsum。只计算一次 1/Wsum，三个方向复用。
    const double inverseWeightSum = 1.0 / weightSum;
    for (int a = 0; a < numElementNodes(); ++a) {
        result.values(a) = inverseWeightSum * weighted(a); // & R = W / S_{W}
        // & R_A = S_{W}^{-1} * (A - R * S_{A})
        // & R_B = S_{W}^{-1} * (B - R * S_{B})
        // & R_C = S_{W}^{-1} * (C - R * S_{C})
        result.parametricGradients.row(a) = inverseWeightSum
            * (weightedDerivatives.row(a) - result.values(a) * derivativeSums.transpose()); // & R_A = S_{W}^{-1} * (A - R * S_{A})
    }
    return result;
}

PointShapeFunctionData TensorProductNurbsSpace::evaluateShapeFunctionsAt(
    const std::array<double, 3>& parameters) const
{
    // 后处理点不一定是 Gauss 点，因此这里才使用 FindSpan 并现场计算一维基函数。
    std::array<int, 3> pointSpans{};
    std::array<Eigen::MatrixXd, 3> derivatives;
    for (int direction = 0; direction < 3; ++direction) {
        pointSpans[direction] = FindSpan(controlPointCounts_[direction] - 1,
                                         degrees_[direction], parameters[direction],
                                         knotVectors_[direction]);
        // &输出矩阵 derivatives 的大小为 (nDer+1) x (p+1)：行表示导数阶数，列表示局部基函数。
        DersBasisFuns(pointSpans[direction], parameters[direction],
                      degrees_[direction], 1, knotVectors_[direction],
                      derivatives[direction]);
    }

    const int d0 = degrees_[0] + 1;
    const int d1 = degrees_[1] + 1;
    const int d2 = degrees_[2] + 1;
    const int base = (pointSpans[0] - degrees_[0])
                   + (pointSpans[1] - degrees_[1]) * controlPointCounts_[0]
                   + (pointSpans[2] - degrees_[2])
                     * controlPointCounts_[0] * controlPointCounts_[1];

    // 与积分点路径相同：先形成 W_A 及其参数导数，再统一执行有理化。
    // values 大小为 n_en，parametricGradients 大小为 n_en x 3。
    PointShapeFunctionData result;
    result.connectivity.resize(static_cast<std::size_t>(numElementNodes()));
    result.shape.values.resize(numElementNodes());
    result.shape.parametricGradients.resize(numElementNodes(), 3);
    Eigen::VectorXd weighted(numElementNodes());
    Eigen::Matrix<double, Eigen::Dynamic, 3> weightedDerivatives(numElementNodes(), 3);
    double weightSum = 0.0;
    Eigen::Vector3d derivativeSums = Eigen::Vector3d::Zero();

    int local = 0;
    for (int a2 = 0; a2 < d2; ++a2) {
        for (int a1 = 0; a1 < d1; ++a1) {
            for (int a0 = 0; a0 < d0; ++a0) {
                const int global = base + a0 + a1 * controlPointCounts_[0]
                                   + a2 * controlPointCounts_[0] * controlPointCounts_[1];
                result.connectivity[static_cast<std::size_t>(local)] = global;
                const double weight = weights_(global);
                const double n0 = derivatives[0](0, a0);
                const double n1 = derivatives[1](0, a1);
                const double n2 = derivatives[2](0, a2);
                weighted(local) = weight * n0 * n1 * n2;
                weightedDerivatives(local, 0) = weight * derivatives[0](1, a0) * n1 * n2;
                weightedDerivatives(local, 1) = weight * n0 * derivatives[1](1, a1) * n2;
                weightedDerivatives(local, 2) = weight * n0 * n1 * derivatives[2](1, a2);
                weightSum += weighted(local);
                derivativeSums += weightedDerivatives.row(local).transpose();
                ++local;
            }
        }
    }

    const double inverseWeightSum = 1.0 / weightSum;
    for (int a = 0; a < numElementNodes(); ++a) {
        result.shape.values(a) = inverseWeightSum * weighted(a);
        result.shape.parametricGradients.row(a) = inverseWeightSum
            * (weightedDerivatives.row(a)
               - result.shape.values(a) * derivativeSums.transpose());
    }
    return result;
}

double TensorProductNurbsSpace::parentJacobianDeterminant(int element) const
{
    const std::array<int, 3> elementSpan = elementSpans(element);
    // 母单元 [-1,1]^3 到参数跨度的映射为 3x3 对角矩阵，
    // 三个对角元分别是各方向跨度长度的一半。
    double determinant = 1.0;
    for (int direction = 0; direction < 3; ++direction) {
        determinant *= (knotVectors_[direction][elementSpan[direction] + 1]
                        - knotVectors_[direction][elementSpan[direction]]) * 0.5;
    }
    return determinant;
}

double TensorProductNurbsSpace::gaussWeight(const std::array<int, 3>& gaussIndices) const
{
    double result = 1.0;
    for (int direction = 0; direction < 3; ++direction) {
        requireIndex(gaussIndices[direction], numGaussPoints(direction),
                     "TensorProductNurbsSpace: Gauss index is out of range");
        result *= gaussRules_[direction].weights[gaussIndices[direction]];
    }
    return result;
}

double TensorProductNurbsSpace::UnivariateTable::value(
    int element, int gauss, int derivative, int basis) const
{
    // 扁平存储布局：[单元][Gauss点][导数阶(0/1)][局部基函数]。
    const std::size_t index = static_cast<std::size_t>(
        (((element * numGauss + gauss) * 2 + derivative) * numBasis + basis));
    return data[index];
}

TensorProductNurbsSpace::GaussRule
TensorProductNurbsSpace::makeGaussLegendreRule(int order)
{
    if (order <= 0) {
        throw std::invalid_argument("Gauss rule order must be positive");
    }
    GaussRule rule;
    rule.points.resize(static_cast<std::size_t>(order));
    rule.weights.resize(static_cast<std::size_t>(order));
    const double pi = std::acos(-1.0);
    // Legendre 多项式的根关于 0 对称，只求一半的根，再镜像到另一侧。
    const int half = (order + 1) / 2;
    for (int i = 0; i < half; ++i) {
        double x = std::cos(pi * (i + 0.75) / (order + 0.5));
        double derivative = 0.0;
        // 用 Newton 法求 P_order(x)=0；current/previous 分别保存相邻两阶多项式。
        for (int iteration = 0; iteration < 100; ++iteration) {
            double previous = 1.0;
            double current = x;
            for (int degree = 2; degree <= order; ++degree) {
                const double next = ((2.0 * degree - 1.0) * x * current
                                     - (degree - 1.0) * previous) / degree;
                previous = current;
                current = next;
            }
            derivative = order * (x * current - previous) / (x * x - 1.0);
            const double step = current / derivative;
            x -= step;
            if (std::abs(step) <= 4.0 * std::numeric_limits<double>::epsilon()) {
                break;
            }
        }
        const double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
        rule.points[static_cast<std::size_t>(i)] = -x;
        rule.points[static_cast<std::size_t>(order - 1 - i)] = x;
        rule.weights[static_cast<std::size_t>(i)] = weight;
        rule.weights[static_cast<std::size_t>(order - 1 - i)] = weight;
    }
    return rule;
}

void TensorProductNurbsSpace::buildSpans()
{
    // 重复结点可能产生零长度跨度；它们不是物理单元，必须跳过。
    for (int direction = 0; direction < 3; ++direction) {
        const int lastControlPoint = controlPointCounts_[direction] - 1;
        for (int i = degrees_[direction]; i <= lastControlPoint; ++i) {
            if (knotVectors_[direction][i] < knotVectors_[direction][i + 1]) {
                spans_[direction].push_back(i);
            }
        }
        if (spans_[direction].empty()) {
            throw std::invalid_argument("TensorProductNurbsSpace: direction has no nonzero knot spans");
        }
    }
}

void TensorProductNurbsSpace::buildTables()
{
    // 每个方向的表布局为 n_elem x n_gauss x 2 x (p+1)，
    // “2”分别对应基函数值和一阶导数。
    for (int direction = 0; direction < 3; ++direction) {
        const int order = degrees_[direction] + 1;
        gaussRules_[direction] = makeGaussLegendreRule(order); // &生成高斯点
        UnivariateTable& table = tables_[direction];
        table.numElements = numElements(direction);
        table.numGauss = order;
        table.numBasis = order;
        table.data.resize(static_cast<std::size_t>(table.numElements * table.numGauss
                                                   * 2 * table.numBasis));
        Eigen::MatrixXd derivatives;
        for (int element = 0; element < table.numElements; ++element) {
            const int knotSpan = spans_[direction][element];
            const double lower = knotVectors_[direction][knotSpan];
            const double upper = knotVectors_[direction][knotSpan + 1];
            for (int gauss = 0; gauss < table.numGauss; ++gauss) {
                // 将 [-1,1] 上的 Gauss 点仿射映射到当前参数跨度 [lower,upper]。 &将母单元的高斯点映射到参数空间
                const double parameter = 0.5 * ((upper - lower) * gaussRules_[direction].points[gauss]
                                                + upper + lower);
                DersBasisFuns(knotSpan, parameter, degrees_[direction], 1,
                              knotVectors_[direction], derivatives);
                for (int derivative = 0; derivative < 2; ++derivative) {
                    for (int basis = 0; basis < table.numBasis; ++basis) {
                        const std::size_t index = static_cast<std::size_t>(
                            (((element * table.numGauss + gauss) * 2 + derivative)
                             * table.numBasis + basis));
                        table.data[index] = derivatives(derivative, basis);
                    }
                }
            }
        }
    }
}

void TensorProductNurbsSpace::validate() const
{
    long long expectedControlPoints = 1;
    for (int direction = 0; direction < 3; ++direction) {
        if (degrees_[direction] < 0 || controlPointCounts_[direction] <= 0
            || !isKnotControlConsistent(knotVectors_[direction].size(),
                                        degrees_[direction], controlPointCounts_[direction])) {
            throw std::invalid_argument("TensorProductNurbsSpace: inconsistent degree and knot vector");
        }
        expectedControlPoints *= controlPointCounts_[direction];
    }
    if (controlPoints_.rows() != expectedControlPoints || controlPoints_.cols() != 3
        || weights_.size() != expectedControlPoints) {
        throw std::invalid_argument("TensorProductNurbsSpace: control net dimensions are inconsistent");
    }
    if (!controlPoints_.allFinite() || !weights_.allFinite()
        || (weights_.array() <= 0.0).any()) {
        throw std::invalid_argument("TensorProductNurbsSpace: coordinates must be finite and weights positive");
    }
}

} // namespace nurbs
