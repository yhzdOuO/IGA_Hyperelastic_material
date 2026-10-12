#ifndef VISUALIZATION_BOUNDARYSAMPLINGCACHE_H
#define VISUALIZATION_BOUNDARYSAMPLINGCACHE_H

#include "visualization/BoundarySampler.h"

namespace visualization {

using SampleScalars = std::array<double, static_cast<int>(ScalarField::Count)>;

/**
 * @brief Precomputes the boundary sampling operators Phi and dPhi/dxi once.
 *
 * Runtime frame updates only multiply the cached basis rows by control-point
 * displacement values; FindSpan and NURBS basis evaluation never enter playback.
 */
class BoundarySamplingCache {
public:
    BoundarySamplingCache(const nurbs::TensorProductNurbsSpace& space,
                          int samplesPerSpan = 4);

    const SurfaceMesh& referenceMesh() const { return referenceMesh_; }

    std::vector<SampleScalars> evaluateScalars(
        const Eigen::VectorXd& displacements,
        const iga::HyperelasticMaterial& material,
        const Eigen::Vector3d& axialDirection = Eigen::Vector3d::UnitX(),
        const Eigen::Vector3d& gravityDirection = -Eigen::Vector3d::UnitZ()) const;

    SurfaceMesh makeFrame(const Eigen::VectorXd& displacements,
                          const std::vector<SampleScalars>& scalars) const;

private:
    struct SampleOperator {
        std::vector<int> connectivity;
        Eigen::VectorXd values;
        iga::GradientMatrix parametricGradients;
        iga::GradientMatrix physicalGradients;
    };

    int dofCount_ = 0;
    SurfaceMesh referenceMesh_;
    std::vector<SampleOperator> operators_;
};

} // namespace visualization

#endif // VISUALIZATION_BOUNDARYSAMPLINGCACHE_H
