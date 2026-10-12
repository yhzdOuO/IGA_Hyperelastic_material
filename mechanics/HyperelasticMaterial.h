#ifndef HYPERELASTICMATERIAL_H
#define HYPERELASTICMATERIAL_H

#include <Eigen/Core>

namespace iga {

using ConstitutiveMatrix = Eigen::Matrix<double, 6, 6>;

/** @brief 一个材料点上完全由当前位移决定的 Total-Lagrangian 运动学状态。 */
struct KinematicState {
    /// 位移梯度 H=du/dX。
    Eigen::Matrix3d displacementGradient = Eigen::Matrix3d::Zero();
    /// 变形梯度 F=I+H。
    Eigen::Matrix3d deformationGradient = Eigen::Matrix3d::Identity();
    /// Green-Lagrange 应变 E=0.5*(F^T*F-I)。
    Eigen::Matrix3d greenLagrangeStrain = Eigen::Matrix3d::Zero();
};

/** @brief 超弹性材料在一个材料点上的应力、相容切线和应变能响应。 */
struct MaterialResponse {
    /// 第二 Piola-Kirchhoff 应力 S。
    Eigen::Matrix3d secondPiolaStress = Eigen::Matrix3d::Zero();
    /**
     * @brief 材料相容切线 D=dS_voigt/dE_engineering。
     *
     * 应力顺序固定为 [S11,S22,S33,S12,S23,S13]，应变顺序固定为
     * [E11,E22,E33,2E12,2E23,2E13]，因此剪切列已经包含工程剪应变因子。
     */
    ConstitutiveMatrix materialTangent = ConstitutiveMatrix::Zero();
    /// 单位参考体积上的应变能密度 Psi。
    double strainEnergyDensity = 0.0;
};

/**
 * @brief 无历史变量的 Total-Lagrangian 超弹性材料通用接口。
 *
 * 求解器只依赖本接口，不依赖具体材料参数。派生类必须在同一次 evaluate()
 * 调用中返回彼此一致的 S、D 和 Psi，确保 Newton 切线与内力来自同一构形。
 */
class HyperelasticMaterial {
public:
    virtual ~HyperelasticMaterial() = default;

    /**
     * @brief 在给定 H、F、E 下计算完整材料响应。
     * @param kinematics 当前材料点的 Total-Lagrangian 运动学量。
     * @return 对应的第二 PK 应力、材料相容切线和应变能密度。
     */
    virtual MaterialResponse evaluate(const KinematicState& kinematics) const = 0;
};

} // namespace iga

#endif // HYPERELASTICMATERIAL_H
