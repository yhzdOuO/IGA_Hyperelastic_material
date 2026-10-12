#ifndef VISUALIZATION_BOUNDARYSAMPLER_H
#define VISUALIZATION_BOUNDARYSAMPLER_H

#include "visualization/FieldEvaluator.h"

namespace visualization {

/** @brief 按六个参数边界面和非零 knot span 生成可渲染三角网格。 */
class BoundarySampler {
public:
    /** @param samplesPerSpan 每个边界单元在两个自由方向上的采样点数。 */
    explicit BoundarySampler(int samplesPerSpan = 4);
    /** @brief 采样六个外表面，同时生成三角形索引和单元边界线索引。 */
    SurfaceMesh sample(const nurbs::TensorProductNurbsSpace& space,
                       const FieldEvaluator& evaluator) const;

private:
    int samplesPerSpan_;
};

} // namespace visualization

#endif // VISUALIZATION_BOUNDARYSAMPLER_H
