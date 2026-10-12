#include "nurbs/NurbsAlgorithm.h"

#include <algorithm>
#include <stdexcept>

namespace nurbs {

// ==================== A2.1：定位结点区间 ================================
int FindSpan(int n, int p, double u, const KnotVector& U)
{
    // 参数必须位于有效参数域 [U[p], U[n+1]]，提前检查可避免二分查找越界。
    if (u < U[p] || u > U[n + 1]) {
        throw std::domain_error("FindSpan: u is outside [U[p], U[n+1]]");
    }

    // 右端点不属于任何半开区间，按 The NURBS Book 的约定归入最后一个跨度。
    if (u >= U[n + 1]) {
        return n;
    }

    // 在可能包含非零 p 次基函数的跨度范围内做二分查找。
    int low  = p;
    int high = n + 1;
    int mid  = (low + high) / 2;
    while (u < U[mid] || u >= U[mid + 1]) {
        if (u < U[mid]) {
            high = mid;
        } else {
            low = mid;
        }
        mid = (low + high) / 2;
    }
    return mid;
}

// ==================== A2.2：非零基函数值 ================================
void BasisFuns(int i, double u, int p, const KnotVector& U,
               std::vector<double>& N)
{
    N.assign(static_cast<std::size_t>(p) + 1, 0.0);

    std::vector<double> left(static_cast<std::size_t>(p) + 1, 0.0);
    std::vector<double> right(static_cast<std::size_t>(p) + 1, 0.0);

    // 逐阶递推：第 j 轮把 j-1 次基函数提升为 j 次基函数。
    // left/right 保存参数点到左右结点的距离，长度均为 p+1。
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j]  = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double temp = N[r] / (right[r + 1] + left[j - r]);
            N[r]  = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        N[j] = saved;
    }
}

// ==================== A2.3：基函数及其导数 ==============================
void DersBasisFuns(int i, double u, int p, int nDer, const KnotVector& U,
                   Eigen::MatrixXd& ders)
{
    if (p < 0 || nDer < 0) {
        throw std::invalid_argument("DersBasisFuns: p and nDer must be nonnegative");
    }
    ders.setZero(nDer + 1, p + 1);

    // 输出矩阵 ders 的大小为 (nDer+1) x (p+1)：行表示导数阶数，列表示局部基函数。
    // p 次多项式的高于 p 阶导数恒为 0，因此实际只计算到 min(nDer,p)。
    const int nd = std::min(nDer, p);

    // ndu 是 (p+1) x (p+1) 的三角递推表，同时保存基函数值和结点差分母。
    std::vector<std::vector<double>> ndu(
        static_cast<std::size_t>(p) + 1,
        std::vector<double>(static_cast<std::size_t>(p) + 1, 0.0));
    std::vector<double> left(static_cast<std::size_t>(p) + 1, 0.0);
    std::vector<double> right(static_cast<std::size_t>(p) + 1, 0.0);

    ndu[0][0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j]  = u - U[i + 1 - j];
        right[j] = U[i + j] - u;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            ndu[j][r] = right[r + 1] + left[j - r];
            const double temp = ndu[r][j - 1] / ndu[j][r];
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    // ndu 最后一列就是 p 次基函数值，写入 ders 的第 0 行。
    for (int j = 0; j <= p; ++j) {
        ders(0, j) = ndu[j][p];
    }

    // a 只有 2 x (p+1)，两行交替使用，避免为每一阶导数重复分配临时数组。
    // 以下对应书中公式 (2.9)，先得到尚未乘阶乘系数的导数。
    std::vector<std::vector<double>> a(
        2, std::vector<double>(static_cast<std::size_t>(p) + 1, 0.0));

    for (int r = 0; r <= p; ++r) {
        int s1 = 0, s2 = 1;
        a[0][0] = 1.0;
        for (int k = 1; k <= nd; ++k) {
            double d = 0.0;
            const int rk = r - k;
            const int pk = p - k;

            if (r >= k) {
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
                d = a[s2][0] * ndu[rk][pk];
            }

            const int j1 = (rk >= -1) ? 1 : -rk;
            const int j2 = (r - 1 <= pk) ? (k - 1) : (p - r);

            for (int j = j1; j <= j2; ++j) {
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
                d += a[s2][j] * ndu[rk + j][pk];
            }

            if (r <= pk) {
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
                d += a[s2][k] * ndu[r][pk];
            }

            ders(k, r) = d;
            std::swap(s1, s2);
        }
    }

    // 最后乘上阶乘系数 p!/(p-k)!，得到真实的 k 阶导数。
    int factor = p;
    for (int k = 1; k <= nd; ++k) {
        for (int j = 0; j <= p; ++j) {
            ders(k, j) *= factor;
        }
        factor *= (p - k);
    }
}

} // namespace nurbs
