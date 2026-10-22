#ifndef NEOHOOKEANMATERIAL_H
#define NEOHOOKEANMATERIAL_H

#include "mechanics/HyperelasticMaterial.h"

namespace iga {

/** @brief 可压缩 Neo-Hookean 材料，使用空间构形公式计算响应。 */
class NeoHookeanMaterial final : public HyperelasticMaterial {
public:
    /** @brief 构造材料，要求 mu>0 且 3*lambda+2*mu>0，与 SVK 参数约定一致。 */
    NeoHookeanMaterial(double lambda, double mu);

    double lambda() const { return lambda_; }
    double mu() const { return mu_; }

    /** @brief 拉回空间响应，得到第二 PK 应力和材料相容切线。 */
    MaterialResponse evaluate(const KinematicState& kinematics) const override;

    /**
     * @brief 由 F 直接计算 Cauchy 应力、空间切线和单位参考体积应变能。
     *
     * Psi=mu/2*(tr(b)-3)-mu*ln(J)+lambda/2*ln(J)^2，b=F*F^T。
     * 空间切线采用 lambda'=lambda/J、mu'=(mu-lambda*ln(J))/J，
     * Voigt 顺序为 [11,22,33,12,23,13]，剪切列采用工程量。
     * 要求 F 有限且 J=det(F)>0。
     */
    SpatialMaterialResponse evaluateSpatial(const KinematicState& kinematics) const override;

private:
    double lambda_;
    double mu_;
};

} // namespace iga

#endif // NEOHOOKEANMATERIAL_H
