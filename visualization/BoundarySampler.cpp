#include "visualization/BoundarySampler.h"

#include <Eigen/Geometry>
#include <array>
#include <stdexcept>

namespace visualization {

BoundarySampler::BoundarySampler(int samplesPerSpan)
    : samplesPerSpan_(samplesPerSpan)
{
    if (samplesPerSpan_ < 2) {
        throw std::invalid_argument("BoundarySampler: samples per span must be at least two");
    }
}

SurfaceMesh BoundarySampler::sample(const nurbs::TensorProductNurbsSpace& space,
                                    const FieldEvaluator& evaluator) const
{
    SurfaceMesh mesh;
    // 功能块 1：每个固定参数方向各有最小、最大两个边界面，共六个面。
    for (int fixedDirection = 0; fixedDirection < 3; ++fixedDirection) {
        // 采用循环顺序指定两个自由参数方向，所有面的局部网格生成逻辑可复用。
        const int free0 = (fixedDirection + 1) % 3;
        const int free1 = (fixedDirection + 2) % 3;
        for (int side = 0; side < 2; ++side) {
            const bool maximumSide = side == 1;
            const double fixedValue = maximumSide
                ? space.knotVector(fixedDirection).backValue()
                : space.knotVector(fixedDirection).frontValue();

            // 功能块 2：在面中心运行时判断 cross(t0,t1) 是否指向外侧。
            // 这避免对六个面手工硬编码翻转表，也适用于一般曲面 patch。
            std::array<double, 3> center{};
            for (int direction = 0; direction < 3; ++direction) {
                center[direction] = 0.5 * (space.knotVector(direction).frontValue()
                                           + space.knotVector(direction).backValue());
            }
            center[fixedDirection] = fixedValue;
            const FieldPoint centerValue = evaluator.evaluate(center);
            const double orientation = centerValue.referenceTangents[free0]
                .cross(centerValue.referenceTangents[free1])
                .dot(centerValue.referenceTangents[fixedDirection]);
            if (std::abs(orientation) <= 1.0e-14) {
                throw std::domain_error("BoundarySampler: degenerate face orientation");
            }
            const bool flip = (maximumSide && orientation < 0.0)
                              || (!maximumSide && orientation > 0.0);
            const double normalSign = flip ? -1.0 : 1.0;

            // 功能块 3：逐非零 knot span 组合采样，不让单个网格格子跨越 knot 线。
            for (int element1 = 0; element1 < space.numElements(free1); ++element1) {
                const int span1 = space.span(free1, element1);
                const double lower1 = space.knotVector(free1)[span1];
                const double upper1 = space.knotVector(free1)[span1 + 1];
                for (int element0 = 0; element0 < space.numElements(free0); ++element0) {
                    const int span0 = space.span(free0, element0);
                    const double lower0 = space.knotVector(free0)[span0];
                    const double upper0 = space.knotVector(free0)[span0 + 1];
                    const std::uint32_t base = static_cast<std::uint32_t>(mesh.vertices.size());

                    // 相邻 span 允许在公共边界生成重复顶点，以保持实现清晰并保留折角。
                    for (int j = 0; j < samplesPerSpan_; ++j) {
                        for (int i = 0; i < samplesPerSpan_; ++i) {
                            std::array<double, 3> parameters = center;
                            parameters[free0] = lower0 + (upper0 - lower0)
                                * static_cast<double>(i) / (samplesPerSpan_ - 1);
                            parameters[free1] = lower1 + (upper1 - lower1)
                                * static_cast<double>(j) / (samplesPerSpan_ - 1);
                            const FieldPoint value = evaluator.evaluate(parameters);
                            // 保存完整 CPU 求值结果；GPU 顶点稍后由 updateGeometry 填充。
                            SurfaceVertex vertex;
                            vertex.parameters = parameters;
                            vertex.referencePosition = value.referencePosition;
                            vertex.displacement = value.displacement;
                            vertex.referenceTangents = value.referenceTangents;
                            vertex.displacementTangents = value.displacementTangents;
                            vertex.scalars = value.scalars;
                            vertex.tangentDirection0 = free0;
                            vertex.tangentDirection1 = free1;
                            vertex.normalSign = normalSign;
                            mesh.vertices.push_back(std::move(vertex));
                        }
                    }
                    // 将 span 内二维采样网格坐标映射为全局顶点索引。
                    const auto index = [this, base](int i, int j) {
                        return base + static_cast<std::uint32_t>(j * samplesPerSpan_ + i);
                    };
                    // 功能块 4：每个四边格使用一致对角线拆为两个同绕序三角形。
                    for (int j = 0; j + 1 < samplesPerSpan_; ++j) {
                        for (int i = 0; i + 1 < samplesPerSpan_; ++i) {
                            const std::uint32_t v00 = index(i, j);
                            const std::uint32_t v10 = index(i + 1, j);
                            const std::uint32_t v11 = index(i + 1, j + 1);
                            const std::uint32_t v01 = index(i, j + 1);
                            // 内向面同时翻转两个三角形的绕序，与解析外法向保持一致。
                            if (!flip) {
                                mesh.triangleIndices.insert(mesh.triangleIndices.end(),
                                                            {v00, v10, v11, v00, v11, v01});
                            } else {
                                mesh.triangleIndices.insert(mesh.triangleIndices.end(),
                                                            {v00, v11, v10, v00, v01, v11});
                            }
                        }
                    }
                    // 功能块 5：线框只记录当前边界单元的四条边，不显示三角形对角线。
                    for (int i = 0; i + 1 < samplesPerSpan_; ++i) {
                        mesh.wireIndices.insert(mesh.wireIndices.end(),
                                                {index(i, 0), index(i + 1, 0),
                                                 index(i, samplesPerSpan_ - 1),
                                                 index(i + 1, samplesPerSpan_ - 1)});
                    }
                    for (int j = 0; j + 1 < samplesPerSpan_; ++j) {
                        mesh.wireIndices.insert(mesh.wireIndices.end(),
                                                {index(0, j), index(0, j + 1),
                                                 index(samplesPerSpan_ - 1, j),
                                                 index(samplesPerSpan_ - 1, j + 1)});
                    }
                }
            }
        }
    }
    // 初始化真实形变倍率和默认轴向应力字段，同时完成 GPU 数据布局和值域统计。
    mesh.updateGeometry(1.0);
    mesh.selectField(ScalarField::AxialStress);
    return mesh;
}

} // namespace visualization
