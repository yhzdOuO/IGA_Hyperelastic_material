#include "mechanics/HyperelasticMaterial.h"

#include <Eigen/LU>
#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

ConstitutiveMatrix voigtStressTransformation(const Eigen::Matrix3d& mapping)
{
    // T 映射应力向量 S -> F*S*F^T。其功共轭映射为 E_engineering=T^T*d_engineering，
    // 所以空间切线是 T*C*T^T/J，不需要额外的剪切倍数。
    constexpr std::array<std::array<int, 2>, 6> pairs{{
        {{0, 0}}, {{1, 1}}, {{2, 2}}, {{0, 1}}, {{1, 2}}, {{0, 2}}
    }};
    ConstitutiveMatrix transformation;
    for (int column = 0; column < 6; ++column) {
        const int i = pairs[column][0];
        const int j = pairs[column][1];
        Eigen::Matrix3d pushed = mapping.col(i) * mapping.col(j).transpose();
        if (i != j) {
            pushed += mapping.col(j) * mapping.col(i).transpose();
        }
        for (int row = 0; row < 6; ++row) {
            transformation(row, column) = pushed(pairs[row][0], pairs[row][1]);
        }
    }
    return transformation;
}

SpatialMaterialResponse pushForwardMaterialResponse(
    const MaterialResponse& response, const Eigen::Matrix3d& deformationGradient)
{
    const double jacobian = deformationGradient.determinant();
    if (!deformationGradient.allFinite() || !std::isfinite(jacobian) || jacobian <= 0.0) {
        throw std::domain_error("pushForwardMaterialResponse: require finite F and det(F)>0");
    }
    if (!response.secondPiolaStress.allFinite() || !response.materialTangent.allFinite()
        || !std::isfinite(response.strainEnergyDensity)) {
        throw std::domain_error("pushForwardMaterialResponse: material response must be finite");
    }

    const ConstitutiveMatrix transformation = voigtStressTransformation(deformationGradient);

    SpatialMaterialResponse result;
    result.cauchyStress = deformationGradient * response.secondPiolaStress
                         * deformationGradient.transpose() / jacobian;
    result.spatialTangent = transformation * response.materialTangent
                            * transformation.transpose() / jacobian;
    result.strainEnergyDensity = response.strainEnergyDensity;
    return result;
}

MaterialResponse pullBackMaterialResponse(
    const SpatialMaterialResponse& response, const Eigen::Matrix3d& deformationGradient)
{
    const double jacobian = deformationGradient.determinant();
    if (!deformationGradient.allFinite() || !std::isfinite(jacobian) || jacobian <= 0.0) {
        throw std::domain_error("pullBackMaterialResponse: require finite F and det(F)>0");
    }
    if (!response.cauchyStress.allFinite() || !response.spatialTangent.allFinite()
        || !std::isfinite(response.strainEnergyDensity)) {
        throw std::domain_error("pullBackMaterialResponse: spatial response must be finite");
    }

    const Eigen::Matrix3d inverse = deformationGradient.inverse();
    const ConstitutiveMatrix transformation = voigtStressTransformation(inverse);
    MaterialResponse result;
    result.secondPiolaStress = jacobian * inverse * response.cauchyStress * inverse.transpose();
    result.materialTangent = jacobian * transformation * response.spatialTangent
                             * transformation.transpose();
    result.strainEnergyDensity = response.strainEnergyDensity;
    return result;
}

} // namespace iga
