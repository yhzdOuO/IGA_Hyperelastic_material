#include "mechanics/NeoHookeanMaterial.h"

#include <Eigen/LU>
#include <cmath>
#include <stdexcept>

namespace iga {

NeoHookeanMaterial::NeoHookeanMaterial(double lambda, double mu)
    : lambda_(lambda), mu_(mu)
{
    if (!std::isfinite(lambda_) || !std::isfinite(mu_) || mu_ <= 0.0
        || 3.0 * lambda_ + 2.0 * mu_ <= 0.0) {
        throw std::invalid_argument(
            "NeoHookeanMaterial: require mu > 0 and 3*lambda + 2*mu > 0");
    }
}

MaterialResponse NeoHookeanMaterial::evaluate(const KinematicState& kinematics) const
{
    return pullBackMaterialResponse(evaluateSpatial(kinematics), kinematics.deformationGradient);
}

SpatialMaterialResponse NeoHookeanMaterial::evaluateSpatial(const KinematicState& kinematics) const
{
    const Eigen::Matrix3d& f = kinematics.deformationGradient;
    const double jacobian = f.determinant();
    if (!f.allFinite() || !std::isfinite(jacobian) || jacobian <= 0.0) {
        throw std::domain_error("NeoHookeanMaterial::evaluateSpatial: require finite F and det(F)>0");
    }

    const Eigen::Matrix3d b = f * f.transpose();
    const double logJ = std::log(jacobian);
    const double effectiveLambda = lambda_ / jacobian;
    const double effectiveMu = (mu_ - lambda_ * logJ) / jacobian;

    // 式 (6.29)：sigma=mu/J*(b-I)+lambda*ln(J)/J*I。
    SpatialMaterialResponse response;
    response.cauchyStress = (mu_ * (b - Eigen::Matrix3d::Identity())
                             + lambda_ * logJ * Eigen::Matrix3d::Identity()) / jacobian;

    // 式 (6.40)-(6.41)：与 SVK 的矩阵形式相同，但有效 Lame 系数随 J 更新。
    // 工程剪切约定下，三个剪切对角项为 mu'，而非 2*mu'。
    response.spatialTangent.topLeftCorner<3, 3>().setConstant(effectiveLambda);
    response.spatialTangent.topLeftCorner<3, 3>().diagonal().array() += 2.0 * effectiveMu;
    response.spatialTangent.bottomRightCorner<3, 3>().diagonal().setConstant(effectiveMu);

    // 式 (6.27)：tr(b)=tr(C)，能量保持按单位参考体积计。
    response.strainEnergyDensity = 0.5 * mu_ * (b.trace() - 3.0) - mu_ * logJ
                                   + 0.5 * lambda_ * logJ * logJ;
    return response;
}

} // namespace iga
