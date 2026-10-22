#ifndef HYPERELASTICMATERIAL_H
#define HYPERELASTICMATERIAL_H

#include <Eigen/Core>

namespace iga {

using ConstitutiveMatrix = Eigen::Matrix<double, 6, 6>;

/** @brief TL/UL 共用的、由总位移决定的材料点运动学状态。 */
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

/** @brief 用于当前构形积分的 Cauchy 应力和空间弹性张量。 */
struct SpatialMaterialResponse {
    Eigen::Matrix3d cauchyStress = Eigen::Matrix3d::Zero();
    /**
     * @brief c=J^{-1} push_forward(C) 的工程剪切 Voigt 矩阵。
     *
     * 顺序为 [11,22,33,12,23,13]。它用于 B^T*D*B，需另加初应力刚度；
     * 不是普通 d(sigma)/d(epsilon)，也不包含 Jaumann 率的应力修正。
     */
    ConstitutiveMatrix spatialTangent = ConstitutiveMatrix::Zero();
    /// 仍为单位参考体积的能量 Psi；当前体积能量密度为 Psi/J。
    double strainEnergyDensity = 0.0;
};

/**
 * @brief 构造应力 Voigt 映射 T，使 voigt(A*S*A^T)=T*voigt(S)。
 *
 * 顺序为 [11,22,33,12,23,13]，应力剪切分量不乘 2。
 * A 可为变形梯度 F 或换基旋转 R；工程剪切切线相应变换为 T*D*T^T。
 * 本函数只执行张量映射，体积因子由调用者处理。
 */
ConstitutiveMatrix voigtStressTransformation(const Eigen::Matrix3d& mapping);

/** @brief 一致推前 S 和 C，要求 F 有限且 det(F)>0。 */
SpatialMaterialResponse pushForwardMaterialResponse(
    const MaterialResponse& response, const Eigen::Matrix3d& deformationGradient);

/**
 * @brief 一致拉回 sigma 和 c，要求 F 有限且 det(F)>0。
 *
 * S=J*F^{-1}*sigma*F^{-T}，C=J*pull_back(c)；能量仍按单位参考体积计。
 */
MaterialResponse pullBackMaterialResponse(
    const SpatialMaterialResponse& response, const Eigen::Matrix3d& deformationGradient);

/**
 * @brief 无历史变量的超弹性材料通用接口，支持参考和空间两种响应。
 *
 * 求解器只依赖本接口，不依赖具体材料参数。派生类必须实现两个响应入口，
 * 可用公式实现其中一个，再用 pushForwardMaterialResponse/pullBackMaterialResponse
 * 实现另一个；每个入口的应力、相容切线和 Psi 必须彼此一致。
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

    /**
     * @brief 计算 Cauchy 应力、空间切线和单位参考体积应变能。
     *
     * 可直接使用空间公式，或推前 evaluate() 的结果；须遵守上述空间切线约定。
     */
    virtual SpatialMaterialResponse evaluateSpatial(const KinematicState& kinematics) const = 0;
};

} // namespace iga

#endif // HYPERELASTICMATERIAL_H
