#ifndef SAINTVENANTKIRCHHOFFMATERIAL_H
#define SAINTVENANTKIRCHHOFFMATERIAL_H

#include "mechanics/HyperelasticMaterial.h"

namespace iga {

/** @brief 用 Lame 参数表示的 Saint Venant-Kirchhoff 超弹性材料。 */
class SaintVenantKirchhoffMaterial final : public HyperelasticMaterial {
public:
    /**
     * @brief 构造各向同性 SVK 材料及常量 6 x 6 Voigt 本构矩阵 D。
     * @param lambda 第一 Lame 参数。
     * @param mu 剪切模量，必须为正。
     */
    SaintVenantKirchhoffMaterial(double lambda, double mu);

    /** @brief 返回第一 Lame 参数 lambda。 */
    double lambda() const { return lambda_; }

    /** @brief 返回剪切模量 mu。 */
    double mu() const { return mu_; }

    /** @brief 返回工程剪应变 Voigt 形式的常量 6 x 6 SVK 切线 D。 */
    const ConstitutiveMatrix& constitutiveMatrix() const { return constitutiveMatrix_; }

    /**
     * @brief 计算 S、常量 D 以及 Psi=0.5*S:E。
     * @param kinematics 当前材料点的 H、F 和 Green-Lagrange 应变 E。
     * @return 与当前构形一致的完整 SVK 材料响应。
     */
    MaterialResponse evaluate(const KinematicState& kinematics) const override;

    /** @brief 一致推前 TL 响应，得到 Cauchy 应力和空间切线。 */
    SpatialMaterialResponse evaluateSpatial(const KinematicState& kinematics) const override;

private:
    double lambda_;
    double mu_;
    ConstitutiveMatrix constitutiveMatrix_ = ConstitutiveMatrix::Zero();
};

} // namespace iga

#endif // SAINTVENANTKIRCHHOFFMATERIAL_H
