#ifndef VISUALIZATION_CANTILEVEREXAMPLE_H
#define VISUALIZATION_CANTILEVEREXAMPLE_H

#include "visualization/SurfaceMesh.h"
#include "nurbs/TensorProductNurbsSpace.h"

namespace visualization {

/** @brief Builds the validated 10x1x1 cubic cantilever analysis patch. */
nurbs::TensorProductNurbsSpace makeCantileverSpace();

/** @brief 示例求解的可视化网格及少量诊断量。 */
struct CantileverVisualization {
    SurfaceMesh mesh;
    double tipDisplacement = 0.0;
    int newtonIterations = 0;
};

/**
 * @brief 构造并求解 tests 中同参数的自重悬臂梁，再生成六面表面网格。
 * @param samplesPerSpan 每个边界 knot span、每个自由方向的采样点数。
 */
CantileverVisualization solveCantileverVisualization(int samplesPerSpan = 4);

} // namespace visualization

#endif // VISUALIZATION_CANTILEVEREXAMPLE_H
