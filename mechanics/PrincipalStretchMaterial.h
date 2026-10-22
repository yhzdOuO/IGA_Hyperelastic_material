#ifndef PRINCIPALSTRETCHMATERIAL_H
#define PRINCIPALSTRETCHMATERIAL_H

#include "mechanics/HyperelasticMaterial.h"

namespace iga {

/**
 * @brief 主伸长对数的二次型各向同性材料，见式 (6.92)-(6.95)。
 *
 * Psi=mu*sum(log(lambda_a)^2)+lambda/2*log(J)^2。
 * lambda_a 为主伸长，材料参数 lambda 为第一 Lame 参数。
 */
class PrincipalStretchMaterial final : public HyperelasticMaterial {
public:
    /** @brief 要求参数有限、mu>0 且 3*lambda+2*mu>0。 */
    PrincipalStretchMaterial(double lambda, double mu);

    double lambda() const { return lambda_; }
    double mu() const { return mu_; }

    /** @brief 拉回空间响应，得到第二 PK 应力和材料相容切线。 */
    MaterialResponse evaluate(const KinematicState& kinematics) const override;

    /**
     * @brief 在 b 的主方向基下计算响应，再用旋转的 Voigt 映射转回笛卡尔系。
     *
     * 要求 F 有限且 det(F)>0；重根采用式 (6.91) 的极限。
     * 空间切线采用工程剪切约定，能量按单位参考体积计。
     */
    SpatialMaterialResponse evaluateSpatial(const KinematicState& kinematics) const override;

private:
    double lambda_;
    double mu_;
};

} // namespace iga

#endif // PRINCIPALSTRETCHMATERIAL_H
