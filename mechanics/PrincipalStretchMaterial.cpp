#include "mechanics/PrincipalStretchMaterial.h"

#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <array>
#include <cmath>
#include <stdexcept>

namespace iga {

PrincipalStretchMaterial::PrincipalStretchMaterial(double lambda, double mu)
    : lambda_(lambda), mu_(mu)
{
    if (!std::isfinite(lambda_) || !std::isfinite(mu_) || mu_ <= 0.0
        || 3.0 * lambda_ + 2.0 * mu_ <= 0.0) {
        throw std::invalid_argument(
            "PrincipalStretchMaterial: require mu > 0 and 3*lambda + 2*mu > 0");
    }
}

MaterialResponse PrincipalStretchMaterial::evaluate(const KinematicState& kinematics) const
{
    return pullBackMaterialResponse(evaluateSpatial(kinematics), kinematics.deformationGradient);
}

SpatialMaterialResponse PrincipalStretchMaterial::evaluateSpatial(const KinematicState& kinematics) const
{
    const Eigen::Matrix3d& f = kinematics.deformationGradient;
    const double jacobian = f.determinant();
    if (!f.allFinite() || !std::isfinite(jacobian) || jacobian <= 0.0) {
        throw std::domain_error(
            "PrincipalStretchMaterial::evaluateSpatial: require finite F and det(F)>0");
    }

    const Eigen::Matrix3d b = f * f.transpose();
    if (!b.allFinite()) {
        throw std::domain_error("PrincipalStretchMaterial::evaluateSpatial: b must be finite");
    }
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> spectral(b);
    if (spectral.info() != Eigen::Success || !spectral.eigenvalues().allFinite()
        || spectral.eigenvalues().minCoeff() <= 0.0) {
        throw std::domain_error(
            "PrincipalStretchMaterial::evaluateSpatial: require positive eigenvalues of b");
    }
    const Eigen::Vector3d squaredStretches = spectral.eigenvalues();
    const Eigen::Vector3d logStretches = 0.5 * squaredStretches.array().log().matrix();
    const double logJ = std::log(jacobian);
    const double scaledMu = mu_ / jacobian;
    // 式 (6.94)：主 Cauchy 应力。
    const Eigen::Vector3d principalStress =
        2.0 * scaledMu * logStretches + (lambda_ * logJ / jacobian) * Eigen::Vector3d::Ones();

    // 式 (6.90)、(6.95)：法向块为 Psi 的对数伸长 Hessian/J，减去 2*sigma_a。
    ConstitutiveMatrix principalTangent = ConstitutiveMatrix::Zero();
    principalTangent.topLeftCorner<3, 3>().setConstant(lambda_ / jacobian);
    principalTangent.topLeftCorner<3, 3>().diagonal() +=
        (2.0 * scaledMu * Eigen::Vector3d::Ones() - 2.0 * principalStress);
    constexpr std::array<std::array<int, 2>, 3> shearPairs{{{{0, 1}}, {{1, 2}}, {{0, 2}}}};
    for (int shear = 0; shear < 3; ++shear) {
        const int a = shearPairs[shear][0], beta = shearPairs[shear][1];
        const double difference = squaredStretches(a) - squaredStretches(beta);
        const double sum = squaredStretches(a) + squaredStretches(beta);
        double coefficient;
        if (std::abs(difference) <= 1e-4 * sum) {
            // 式 (6.91)：重根时为 mu/J-sigma_a。
            // 近重根用等价式 (mu/J)*atanh(x)/x-(sigma_a+sigma_b)/2 的级数，
            // x=(lambda_a^2-lambda_b^2)/(lambda_a^2+lambda_b^2)，避免相减消去。
            const double x = difference / sum;
            const double x2 = x * x;
            coefficient = scaledMu * (1.0 + x2 / 3.0 + x2 * x2 / 5.0)
                          - 0.5 * (principalStress(a) + principalStress(beta));
        } else {
            coefficient = (principalStress(a) * squaredStretches(beta)
                           - principalStress(beta) * squaredStretches(a)) / difference;
        }
        principalTangent(3 + shear, 3 + shear) = coefficient;
    }

    // R 的列为当前主方向 n_a；调整符号使其为正旋转，不改变张量响应。
    Eigen::Matrix3d rotation = spectral.eigenvectors();
    if (rotation.determinant() < 0.0) rotation.col(0) *= -1.0;
    const ConstitutiveMatrix transformation = voigtStressTransformation(rotation);
    Eigen::Matrix<double, 6, 1> principalVoigt = Eigen::Matrix<double, 6, 1>::Zero();
    principalVoigt.head<3>() = principalStress;
    const Eigen::Matrix<double, 6, 1> cartesianStress = transformation * principalVoigt;

    SpatialMaterialResponse response;
    response.cauchyStress << cartesianStress(0), cartesianStress(3), cartesianStress(5),
                             cartesianStress(3), cartesianStress(1), cartesianStress(4),
                             cartesianStress(5), cartesianStress(4), cartesianStress(2);
    // 这是旋转换基，主方向切线已经包含 1/J，此处不再乘体积因子。
    response.spatialTangent = transformation * principalTangent * transformation.transpose();
    response.strainEnergyDensity = mu_ * logStretches.squaredNorm() + 0.5 * lambda_ * logJ * logJ;
    return response;
}

} // namespace iga
