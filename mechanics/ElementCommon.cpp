#include "mechanics/ElementCommon.h"

#include <stdexcept>

namespace iga {

BMatrix makeLinearB(const Eigen::Vector3d& physicalGradient)
{
    const double g1 = physicalGradient(0);
    const double g2 = physicalGradient(1);
    const double g3 = physicalGradient(2);
    BMatrix result;
    result << g1, 0.0, 0.0,
              0.0, g2, 0.0,
              0.0, 0.0, g3,
              g2, g1, 0.0,
              0.0, g3, g2,
              g3, 0.0, g1;
    return result;
}

KinematicState evaluateKinematics(
    const nurbs::ControlPointMatrix& localDisplacements,
    const GradientMatrix& referenceGradients)
{
    if (localDisplacements.rows() != referenceGradients.rows()
        || localDisplacements.cols() != 3 || referenceGradients.cols() != 3
        || !localDisplacements.allFinite() || !referenceGradients.allFinite()) {
        throw std::invalid_argument("evaluateKinematics: displacement dimensions are inconsistent");
    }
    KinematicState result;
    result.displacementGradient.noalias() = localDisplacements.transpose()
                                            * referenceGradients;
    result.deformationGradient = Eigen::Matrix3d::Identity()
                                 + result.displacementGradient;
    // 用 H 形式避免小位移下 F^T*F-I 的浮点抵消；TL/UL 使用同一总变形。
    const Eigen::Matrix3d quadratic = result.displacementGradient.transpose()
                                      * result.displacementGradient;
    result.greenLagrangeStrain = 0.5 * (result.displacementGradient
                                        + result.displacementGradient.transpose()
                                        + quadratic);
    return result;
}

Eigen::Matrix<double, 6, 1> stressVoigt(const Eigen::Matrix3d& stress)
{
    Eigen::Matrix<double, 6, 1> result;
    result << stress(0, 0), stress(1, 1), stress(2, 2),
              stress(0, 1), stress(1, 2), stress(0, 2);
    return result;
}

} // namespace iga
