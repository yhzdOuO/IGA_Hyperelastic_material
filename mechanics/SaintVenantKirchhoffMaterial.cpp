#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <cmath>
#include <stdexcept>

namespace iga {

SaintVenantKirchhoffMaterial::SaintVenantKirchhoffMaterial(double lambda, double mu)
    : lambda_(lambda), mu_(mu)
{
    if (!std::isfinite(lambda_) || !std::isfinite(mu_) || mu_ <= 0.0
        || 3.0 * lambda_ + 2.0 * mu_ <= 0.0) {
        throw std::invalid_argument(
            "SaintVenantKirchhoffMaterial: require mu > 0 and 3*lambda + 2*mu > 0");
    }

    // ---------- 功能块 1：预计算与构形无关的 SVK 材料切线 ----------
    // Voigt 顺序固定为 [11,22,33,12,23,13]；剪切应变使用工程量，
    // 所以三个剪切对角项是 mu，而不是张量表达中的 2*mu。
    constitutiveMatrix_.setZero();
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            constitutiveMatrix_(row, column) = lambda_;
        }
        constitutiveMatrix_(row, row) += 2.0 * mu_;
    }
    constitutiveMatrix_(3, 3) = mu_;
    constitutiveMatrix_(4, 4) = mu_;
    constitutiveMatrix_(5, 5) = mu_;
}

MaterialResponse SaintVenantKirchhoffMaterial::evaluate(
    const KinematicState& kinematics) const
{
    if (!kinematics.displacementGradient.allFinite()
        || !kinematics.deformationGradient.allFinite()
        || !kinematics.greenLagrangeStrain.allFinite()) {
        throw std::invalid_argument(
            "SaintVenantKirchhoffMaterial::evaluate: kinematics must be finite");
    }

    // ---------- 功能块 1：由 Green-Lagrange 应变计算第二 PK 应力 ----------
    // E 和 S 都保留 3x3 张量形式，避免在本构公式中混入 Voigt 剪切因子。
    const Eigen::Matrix3d& strain = kinematics.greenLagrangeStrain;
    MaterialResponse response;
    response.secondPiolaStress =
        lambda_ * strain.trace() * Eigen::Matrix3d::Identity()
        + 2.0 * mu_ * strain;

    // ---------- 功能块 2：提交相容切线和单位参考体积应变能 ----------
    // SVK 的 D 与构形无关；Psi=0.5*S:E 中矩阵逐项积会自然计入两次剪切项。
    response.materialTangent = constitutiveMatrix_;
    response.strainEnergyDensity = 0.5
        * response.secondPiolaStress.cwiseProduct(strain).sum();
    return response;
}

} // namespace iga
