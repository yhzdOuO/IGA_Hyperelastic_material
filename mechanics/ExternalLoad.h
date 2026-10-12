#ifndef EXTERNALLOAD_H
#define EXTERNALLOAD_H

#include <Eigen/Core>
#include <functional>

namespace iga {

/** @brief 返回指定物理时刻完整全局外载荷向量的回调类型；空回调表示零载荷。 */
using TimeDependentLoad = std::function<Eigen::VectorXd(double)>;

} // namespace iga

#endif // EXTERNALLOAD_H
