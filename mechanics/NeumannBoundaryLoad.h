#ifndef NEUMANNBOUNDARYLOAD_H
#define NEUMANNBOUNDARYLOAD_H

#include "mechanics/ExternalLoad.h"
#include "nurbs/TensorProductNurbsSpace.h"

#include <functional>

namespace iga {

/**
 * @brief 返回单位参考面积牵引力向量 T0(t) 的时间回调。
 *
 * 返回值在全局参考坐标系中同时表达大小和方向，不表示整个面的合力。
 */
using ReferenceTraction = std::function<Eigen::Vector3d(double)>;

/**
 * @brief 为一个未裁剪参数边界面创建参考构形 Neumann 牵引力载荷。
 *
 * 创建时预积分 c_A=integral_Gamma0(R_A dA0)，随后返回的时间回调只计算
 * f_A(t)=c_A*T0(t)。该载荷是方向固定于参考坐标系的死载荷，不包含 follower
 * load 所需的当前法向、当前面积或额外载荷切线。
 *
 * @param space 三维 tensor-product NURBS patch 及参考控制点。
 * @param normalDirection 固定参数方向；0=xi、1=eta、2=zeta。
 * @param side 参数面的端点；0=最小端、1=最大端。
 * @param traction 返回单位参考面积牵引力 T0(t) 的非空回调。
 * @return 可直接传给动力学求解器的 N_dof 维时间载荷回调。
 */
TimeDependentLoad makeReferenceFaceTractionLoad(
    const nurbs::TensorProductNurbsSpace& space,
    int normalDirection,
    int side,
    ReferenceTraction traction);

} // namespace iga

#endif // NEUMANNBOUNDARYLOAD_H
