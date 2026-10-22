#include "mechanics/UpdatedLagrangianElement.h"
#include "mechanics/TotalLagrangianElement.h"
#include "mechanics/GeneralizedAlphaDynamics.h"
#include "mechanics/SaintVenantKirchhoffMaterial.h"
#include "mechanics/NeoHookeanMaterial.h"
#include "mechanics/PrincipalStretchMaterial.h"

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using nurbs::ControlPointMatrix;
using nurbs::KnotVector;
using nurbs::TensorProductNurbsSpace;
constexpr std::array<std::array<int, 2>, 6> pairs{{
    {{0, 0}}, {{1, 1}}, {{2, 2}}, {{0, 1}}, {{1, 2}}, {{0, 2}}
}};

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

template<class A, class B>
void requireClose(const A& actual, const B& expected, double tolerance,
                  const std::string& message)
{
    const double error = (actual - expected).norm() / std::max(1.0, expected.norm());
    if (error > tolerance || !std::isfinite(error)) {
        throw std::runtime_error(message + ": relative error=" + std::to_string(error));
    }
}

template<class Function>
void requireDomainError(Function function, const std::string& message)
{
    bool threw = false;
    try { function(); }
    catch (const std::domain_error&) { threw = true; }
    require(threw, message);
}

// 非单位参数区间、多 span、非均匀权重和弯曲几何，避免母单元因子被漏乘后仍能通过。
TensorProductNurbsSpace makePatch()
{
    const std::array<std::vector<double>, 3> knots{{
        {-2.0, -2.0, -2.0, 0.0, 3.0, 3.0, 3.0},
        {1.0, 1.0, 1.0, 2.0, 5.0, 5.0, 5.0},
        {-1.0, -1.0, -1.0, 2.0, 2.0, 2.0}
    }};
    ControlPointMatrix points(48, 3);
    Eigen::VectorXd weights(48);
    int index = 0;
    for (int k = 0; k < 3; ++k) {
        for (int j = 0; j < 4; ++j) {
            for (int i = 0; i < 4; ++i) {
                const double x = (0.5 * (knots[0][i + 1] + knots[0][i + 2]) + 2.0) / 5.0;
                const double y = (0.5 * (knots[1][j + 1] + knots[1][j + 2]) - 1.0) / 4.0;
                const double z = (0.5 * (knots[2][k + 1] + knots[2][k + 2]) + 1.0) / 3.0;
                points.row(index) << 1.4*x + 0.08*y*z, 0.9*y + 0.06*x*z, 1.1*z + 0.05*x*y;
                weights(index) = 1.0 + 0.08 * std::sin(0.7 * index);
                ++index;
            }
        }
    }
    return TensorProductNurbsSpace(
        {KnotVector(knots[0]), KnotVector(knots[1]), KnotVector(knots[2])},
        {2, 2, 2}, std::move(points), std::move(weights));
}

Eigen::VectorXd makeDisplacements(const TensorProductNurbsSpace& space)
{
    Eigen::Matrix3d stretch;
    stretch << 1.12, 0.12, -0.04, 0.03, 0.91, 0.08, 0.06, -0.02, 1.05;
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        0.45, Eigen::Vector3d(1.0, 2.0, -1.0).normalized()).toRotationMatrix();
    Eigen::VectorXd result(iga::numDofs(space));
    for (int a = 0; a < space.controlPoints().rows(); ++a) {
        const Eigen::Vector3d x = space.controlPoints().row(a).transpose();
        const Eigen::Vector3d perturbation(0.04*x(1)*x(2), -0.03*x(0)*x(2), 0.05*x(0)*x(1));
        result.segment<3>(3*a) = rotation*(stretch*x + perturbation) - x
                                + Eigen::Vector3d(0.2, -0.1, 0.3);
    }
    return result;
}

// 测试用双重公式：TL 从 C 求解，UL 直接采用书中 9.39 的空间表达。
// 用它验证材料可以覆盖 evaluateSpatial()，且 UL 确实使用该入口。
class DirectSpatialNeoHookean final : public iga::HyperelasticMaterial {
public:
    mutable int spatialCalls = 0;

    iga::MaterialResponse evaluate(const iga::KinematicState& state) const override
    {
        const Eigen::Matrix3d c = state.deformationGradient.transpose() * state.deformationGradient;
        const Eigen::Matrix3d inverse = c.inverse();
        const double logJ = std::log(state.deformationGradient.determinant());
        iga::MaterialResponse result;
        result.secondPiolaStress = mu*(Eigen::Matrix3d::Identity() - inverse) + lambda*logJ*inverse;
        for (int r = 0; r < 6; ++r) {
            const int i = pairs[r][0], j = pairs[r][1];
            for (int s = 0; s < 6; ++s) {
                const int k = pairs[s][0], l = pairs[s][1];
                result.materialTangent(r,s) = lambda*inverse(i,j)*inverse(k,l)
                    + (mu-lambda*logJ)*(inverse(i,k)*inverse(j,l)+inverse(i,l)*inverse(j,k));
            }
        }
        result.strainEnergyDensity = 0.5*mu*(c.trace()-3.0)-mu*logJ+0.5*lambda*logJ*logJ;
        return result;
    }

    iga::SpatialMaterialResponse evaluateSpatial(const iga::KinematicState& state) const override
    {
        ++spatialCalls;
        const Eigen::Matrix3d b = state.deformationGradient * state.deformationGradient.transpose();
        const double jacobian = state.deformationGradient.determinant();
        const double logJ = std::log(jacobian);
        const double lp = lambda/jacobian, mp = (mu-lambda*logJ)/jacobian;
        iga::SpatialMaterialResponse result;
        result.cauchyStress = (mu*(b-Eigen::Matrix3d::Identity())
                              + lambda*logJ*Eigen::Matrix3d::Identity())/jacobian;
        result.spatialTangent.topLeftCorner<3,3>().setConstant(lp);
        result.spatialTangent.topLeftCorner<3,3>().diagonal().array() += 2.0*mp;
        result.spatialTangent.bottomRightCorner<3,3>().diagonal().setConstant(mp);
        result.strainEnergyDensity = 0.5*mu*(b.trace()-3.0)-mu*logJ+0.5*lambda*logJ*logJ;
        return result;
    }
private:
    const double lambda = 5.0, mu = 2.0;
};

void testSpatialPushForward()
{
    Eigen::Matrix3d f;
    f << 1.1, 0.3, -0.1, 0.2, 0.9, 0.15, 0.05, -0.08, 1.2;
    iga::MaterialResponse reference;
    iga::ConstitutiveMatrix m;
    for (int r = 0; r < 6; ++r)
        for (int s = 0; s < 6; ++s) m(r,s) = std::sin(0.9*r+0.7*s) + (r==s ? 2.0 : 0.0);
    reference.materialTangent = m.transpose()*m;
    reference.secondPiolaStress << 2.0, 0.3, -0.2, 0.3, -1.0, 0.7, -0.2, 0.7, 0.5;
    reference.strainEnergyDensity = 1.25;
    const auto actual = iga::pushForwardMaterialResponse(reference, f);
    const auto voigtIndex = [](int i, int j) {
        if (i==j) return i;
        if ((i==0 && j==1)||(i==1 && j==0)) return 3;
        if ((i==1 && j==2)||(i==2 && j==1)) return 4;
        return 5;
    };
    iga::ConstitutiveMatrix expected = iga::ConstitutiveMatrix::Zero();
    for (int r = 0; r < 6; ++r) for (int s = 0; s < 6; ++s)
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
                expected(r,s) += f(pairs[r][0],i)*f(pairs[r][1],j)
                    * reference.materialTangent(voigtIndex(i,j),voigtIndex(k,l))
                    * f(pairs[s][0],k)*f(pairs[s][1],l)/f.determinant();
    requireClose(actual.spatialTangent, expected, 2e-14, "anisotropic fourth-order push-forward");
    requireClose(actual.cauchyStress, (f*reference.secondPiolaStress*f.transpose()/f.determinant()).eval(),
                 2e-14, "Cauchy stress push-forward");
    require(actual.strainEnergyDensity==reference.strainEnergyDensity, "energy remains per reference volume");

    const auto recovered = iga::pullBackMaterialResponse(actual, f);
    requireClose(recovered.secondPiolaStress, reference.secondPiolaStress, 2e-14, "stress round trip");
    requireClose(recovered.materialTangent, reference.materialTangent, 2e-14, "anisotropic tangent round trip");
    require(recovered.strainEnergyDensity==reference.strainEnergyDensity, "pull-back preserves reference energy");
}

void testSpatialPullBack()
{
    iga::KinematicState state;
    state.deformationGradient << 1.1, 0.3, -0.1, 0.2, 0.9, 0.15, 0.05, -0.08, 1.2;
    const DirectSpatialNeoHookean material;
    const auto spatial = material.evaluateSpatial(state);
    const auto expected = material.evaluate(state);
    const auto actual = iga::pullBackMaterialResponse(spatial, state.deformationGradient);
    requireClose(actual.secondPiolaStress, expected.secondPiolaStress, 2e-14, "Neo-Hookean stress pull-back");
    requireClose(actual.materialTangent, expected.materialTangent, 2e-14, "Neo-Hookean tangent pull-back");
    require(std::abs(actual.strainEnergyDensity-expected.strainEnergyDensity)<2e-14, "Neo-Hookean energy");

    Eigen::Matrix3d invalid = Eigen::Matrix3d::Identity();
    invalid(0,0) = 0.0;
    requireDomainError([&] { iga::pullBackMaterialResponse(spatial, invalid); }, "reject singular F");
    invalid(0,0) = -1.0;
    requireDomainError([&] { iga::pullBackMaterialResponse(spatial, invalid); }, "reject inverted F");
    invalid(0,0) = std::numeric_limits<double>::quiet_NaN();
    requireDomainError([&] { iga::pullBackMaterialResponse(spatial, invalid); }, "reject nonfinite F");

    auto invalidResponse = spatial;
    invalidResponse.cauchyStress(0,0) = std::numeric_limits<double>::infinity();
    requireDomainError([&] { iga::pullBackMaterialResponse(invalidResponse, state.deformationGradient); },
                       "reject nonfinite stress");
    invalidResponse = spatial;
    invalidResponse.spatialTangent(3,4) = std::numeric_limits<double>::infinity();
    requireDomainError([&] { iga::pullBackMaterialResponse(invalidResponse, state.deformationGradient); },
                       "reject nonfinite tangent");
    invalidResponse = spatial;
    invalidResponse.strainEnergyDensity = std::numeric_limits<double>::infinity();
    requireDomainError([&] { iga::pullBackMaterialResponse(invalidResponse, state.deformationGradient); },
                       "reject nonfinite energy");
}

void testNeoHookeanMaterial()
{
    const iga::NeoHookeanMaterial material(5.0, 2.0);
    const DirectSpatialNeoHookean reference;
    const iga::SaintVenantKirchhoffMaterial svk(5.0, 2.0);
    const iga::KinematicState initial;
    const auto initialResponse = material.evaluateSpatial(initial);
    require(initialResponse.cauchyStress.isZero() && initialResponse.strainEnergyDensity == 0.0,
            "Neo-Hookean initial stress and energy");
    requireClose(initialResponse.spatialTangent, svk.constitutiveMatrix(), 2e-14,
                 "Neo-Hookean initial tangent equals linear elasticity");

    Eigen::Matrix3d mixed;
    mixed << 1.1, 0.3, -0.1, 0.2, 0.9, 0.15, 0.05, -0.08, 1.2;
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        0.7, Eigen::Vector3d(1.0, -2.0, 3.0).normalized()).toRotationMatrix();
    const std::array<Eigen::Matrix3d, 4> gradients{{
        0.8 * Eigen::Matrix3d::Identity(), 1.3 * Eigen::Matrix3d::Identity(), mixed, rotation
    }};
    for (const auto& f : gradients) {
        iga::KinematicState state;
        state.deformationGradient = f;
        const auto expected = reference.evaluate(state);
        const auto actual = material.evaluate(state);
        const auto spatial = material.evaluateSpatial(state);
        const auto pushed = iga::pushForwardMaterialResponse(expected, f);
        requireClose(actual.secondPiolaStress, expected.secondPiolaStress, 2e-14, "Neo-Hookean TL stress");
        requireClose(actual.materialTangent, expected.materialTangent, 2e-14, "Neo-Hookean TL tangent");
        requireClose(spatial.cauchyStress, pushed.cauchyStress, 2e-14, "Neo-Hookean spatial stress");
        requireClose(spatial.spatialTangent, pushed.spatialTangent, 2e-14, "Neo-Hookean spatial tangent");
        require(std::abs(spatial.strainEnergyDensity-expected.strainEnergyDensity)<2e-14,
                "Neo-Hookean reference-volume energy");

        const double h = 1e-6;
        auto plus = state, minus = state;
        plus.deformationGradient += h * mixed;
        minus.deformationGradient -= h * mixed;
        const double energyDerivative = (material.evaluateSpatial(plus).strainEnergyDensity
                                        - material.evaluateSpatial(minus).strainEnergyDensity) / (2.0*h);
        const double stressPower = (f * actual.secondPiolaStress).cwiseProduct(mixed).sum();
        require(std::abs(energyDerivative-stressPower)<2e-8*std::max(1.0,std::abs(stressPower)),
                "Neo-Hookean stress is the energy derivative");
    }

    for (double value : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()}) {
        iga::KinematicState state;
        state.deformationGradient(0,0) = value;
        requireDomainError([&] { material.evaluateSpatial(state); }, "Neo-Hookean rejects invalid F in UL");
        requireDomainError([&] { material.evaluate(state); }, "Neo-Hookean rejects invalid F in TL");
    }
    const std::array<std::array<double, 2>, 4> invalidParameters{{
        {{5.0, 0.0}}, {{5.0, -1.0}}, {{-2.0, 2.0}},
        {{std::numeric_limits<double>::infinity(), 2.0}}
    }};
    for (const auto& parameters : invalidParameters) {
        bool threw = false;
        try { iga::NeoHookeanMaterial invalid(parameters[0], parameters[1]); }
        catch (const std::invalid_argument&) { threw = true; }
        require(threw, "Neo-Hookean rejects invalid Lame parameters");
    }
}

void testPrincipalStretchMaterial()
{
    const iga::PrincipalStretchMaterial material(5.0, 2.0);
    const iga::SaintVenantKirchhoffMaterial svk(5.0, 2.0);
    const auto initial = material.evaluateSpatial(iga::KinematicState{});
    require(initial.cauchyStress.isZero() && initial.strainEnergyDensity == 0.0,
            "principal-stretch initial stress and energy");
    requireClose(initial.spatialTangent, svk.constitutiveMatrix(), 2e-14,
                 "principal-stretch initial tangent");

    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        0.7, Eigen::Vector3d(1.0, -2.0, 3.0).normalized()).toRotationMatrix();
    const Eigen::Matrix3d rightRotation = Eigen::AngleAxisd(
        -0.4, Eigen::Vector3d(2.0, 1.0, -1.0).normalized()).toRotationMatrix();
    const std::array<Eigen::Vector3d, 7> stretches{{
        {1.0, 1.0, 1.0}, {0.8, 0.8, 0.8}, {1.3, 1.3, 1.3},
        {1.2, 1.2, 0.9}, {1.2, 1.2*(1.0+1e-10), 0.9},
        {1.2, 1.2003, 0.9}, {0.7, 1.1, 1.6}
    }};
    for (const auto& stretch : stretches) {
        iga::KinematicState state;
        state.deformationGradient = rotation * stretch.asDiagonal() * rightRotation.transpose();
        const double jacobian = stretch.prod();
        const double logJ = std::log(jacobian);
        const Eigen::Vector3d logs = stretch.array().log().matrix();
        const Eigen::Vector3d principalStress =
            (4.0 * logs + 5.0 * logJ * Eigen::Vector3d::Ones()) / jacobian;
        const auto spatial = material.evaluateSpatial(state);
        requireClose(spatial.cauchyStress,
                     (rotation * principalStress.asDiagonal() * rotation.transpose()).eval(),
                     3e-14, "principal-stretch stress in Cartesian coordinates");
        const double expectedEnergy = 2.0 * logs.squaredNorm() + 2.5 * logJ * logJ;
        require(std::abs(spatial.strainEnergyDensity-expectedEnergy)<3e-14,
                "principal-stretch reference-volume energy");
        const auto reference = material.evaluate(state);
        const auto pushed = iga::pushForwardMaterialResponse(reference, state.deformationGradient);
        requireClose(pushed.cauchyStress, spatial.cauchyStress, 3e-14, "principal-stretch TL/UL stress");
        requireClose(pushed.spatialTangent, spatial.spatialTangent, 3e-14, "principal-stretch TL/UL tangent");

        // 独立验证完整六列：sigma_dot=c:d+L*sigma+sigma*L^T-tr(L)*sigma。
        // 在旋转后的三重根、双重根、近重根及不同根上均检查工程剪切系数。
        for (int column = 0; column < 6; ++column) {
            Eigen::Matrix3d direction = Eigen::Matrix3d::Zero();
            const int i = pairs[column][0], j = pairs[column][1];
            direction(i,j) = direction(j,i) = (i == j ? 1.0 : 0.5);
            const double h = 1e-6;
            auto plus = state, minus = state;
            plus.deformationGradient = (Eigen::Matrix3d::Identity()+h*direction)*state.deformationGradient;
            minus.deformationGradient = (Eigen::Matrix3d::Identity()-h*direction)*state.deformationGradient;
            const auto plusResponse = material.evaluateSpatial(plus);
            const auto minusResponse = material.evaluateSpatial(minus);
            const Eigen::Matrix3d tangentAction =
                (plusResponse.cauchyStress-minusResponse.cauchyStress)/(2.0*h)
                - direction*spatial.cauchyStress - spatial.cauchyStress*direction.transpose()
                + direction.trace()*spatial.cauchyStress;
            Eigen::Matrix<double, 6, 1> expected;
            for (int row = 0; row < 6; ++row) expected(row) = tangentAction(pairs[row][0],pairs[row][1]);
            requireClose(spatial.spatialTangent.col(column), expected, 2e-8,
                         "principal-stretch spatial tangent finite difference");
            const double energyDerivative =
                (plusResponse.strainEnergyDensity-minusResponse.strainEnergyDensity)/(2.0*h);
            const double stressPower = jacobian*spatial.cauchyStress.cwiseProduct(direction).sum();
            require(std::abs(energyDerivative-stressPower)<2e-8*std::max(1.0,std::abs(stressPower)),
                    "principal-stretch stress is the energy derivative");
        }
    }
}

void testPrincipalStretchInvalidInputs()
{
    const iga::PrincipalStretchMaterial material(5.0, 2.0);
    for (double value : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()}) {
        iga::KinematicState state;
        state.deformationGradient(0,0) = value;
        requireDomainError([&] { material.evaluateSpatial(state); }, "principal-stretch rejects invalid F in UL");
        requireDomainError([&] { material.evaluate(state); }, "principal-stretch rejects invalid F in TL");
    }
    const std::array<std::array<double, 2>, 4> invalidParameters{{
        {{5.0, 0.0}}, {{5.0, -1.0}}, {{-2.0, 2.0}},
        {{std::numeric_limits<double>::infinity(), 2.0}}
    }};
    for (const auto& parameters : invalidParameters) {
        bool threw = false;
        try { iga::PrincipalStretchMaterial invalid(parameters[0], parameters[1]); }
        catch (const std::invalid_argument&) { threw = true; }
        require(threw, "principal-stretch rejects invalid Lame parameters");
    }
}

void testGeometryAndMeasures()
{
    const auto space = makePatch();
    const auto displacement = makeDisplacements(space);
    double minimumJ = 1e10, maximumJ = 0.0;
    for (int e = 0; e < space.numElements(); ++e) {
        const auto points = space.elementControlPoints(e);
        const auto local = iga::gatherElementDisplacements(space,e,displacement);
        const ControlPointMatrix current = points+local;
        for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int i = 0; i < 3; ++i) {
            const std::array<int,3> q{i,j,k};
            const auto shape = space.evaluateShapeFunctions(e,q);
            const auto ref = iga::evaluateReferenceGradients(shape, points,
                space.parentJacobianDeterminant(e)*space.gaussWeight(q));
            const auto point = iga::evaluateSpatialPoint(ref,local,current);
            const auto tl = iga::evaluateReferencePoint(space,e,q);
            const Eigen::Matrix3d& f = point.kinematics.deformationGradient;
            requireClose(ref.referenceGradients,tl.referenceGradients,2e-13,"reference gradients");
            requireClose(point.spatialGradients,(ref.referenceGradients*f.inverse()).eval(),
                         2e-13,"spatial gradient chain rule");
            requireClose(point.currentJacobian,(f*tl.referenceJacobian).eval(),2e-13,"A_t=F*A_0");
            require(std::abs(point.integrationWeight-point.deformationJacobian*tl.integrationWeight)<2e-13,
                    "dv=J*dV with parent-to-parameter weight exactly once");
            require(std::abs(point.deformationJacobian-point.currentJacobianDeterminant
                             /ref.referenceJacobianDeterminant)<2e-13,"J=det(A_t)/det(A_0)");
            minimumJ = std::min(minimumJ,point.deformationJacobian);
            maximumJ = std::max(maximumJ,point.deformationJacobian);
        }
    }
    require(maximumJ-minimumJ>1e-3,"test deformation must be nonuniform across quadrature points");
}

void testElementEquivalenceAndTangent()
{
    const auto space = makePatch();
    const auto displacement = makeDisplacements(space);
    const iga::SaintVenantKirchhoffMaterial svk(5.0,2.0);
    const DirectSpatialNeoHookean neo;
    const iga::NeoHookeanMaterial neoHookean(5.0, 2.0);
    const iga::PrincipalStretchMaterial principalStretch(5.0, 2.0);
    const std::array<const iga::HyperelasticMaterial*,4> materials{&svk,&neo,&neoHookean,&principalStretch};
    for (const auto* material : materials) {
        for (int e = 0; e < space.numElements(); ++e) {
            const auto local = iga::gatherElementDisplacements(space,e,displacement);
            const auto tl = iga::computeTotalLagrangianElementResponse(space,e,local,*material);
            const auto ul = iga::computeUpdatedLagrangianElementResponse(space,e,local,*material);
            requireClose(ul.internalForce,tl.internalForce,3e-12,"TL/UL element internal force");
            requireClose(ul.tangentStiffness,tl.tangentStiffness,3e-12,"TL/UL element tangent");
            requireClose(ul.tangentStiffness,ul.tangentStiffness.transpose().eval(),
                         3e-13,"UL tangent symmetry without forced mirroring");
            ControlPointMatrix direction(local.rows(),3);
            for (int a = 0; a < direction.rows(); ++a)
                for (int i = 0; i < 3; ++i) direction(a,i)=std::sin(0.7*(3*a+i+1));
            direction.normalize();
            Eigen::VectorXd flat(3*direction.rows());
            for (int a = 0; a < direction.rows(); ++a) flat.segment<3>(3*a)=direction.row(a).transpose();
            const double h = 1e-6;
            const auto plus = iga::computeUpdatedLagrangianElementResponse(space,e,(local+h*direction).eval(),*material);
            const auto minus = iga::computeUpdatedLagrangianElementResponse(space,e,(local-h*direction).eval(),*material);
            requireClose(((plus.internalForce-minus.internalForce)/(2*h)).eval(),
                         (ul.tangentStiffness*flat).eval(),3e-8,"UL internal-force finite-difference tangent");
        }
    }
    require(neo.spatialCalls>0,"UL must dispatch to the direct spatial material interface");
}

void testRigidMotionAndInvalidStates()
{
    const auto space = makePatch();
    const iga::SaintVenantKirchhoffMaterial material(5.0,2.0);
    const auto points = space.elementControlPoints(0);
    const Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        1.2,Eigen::Vector3d(1.0,-2.0,3.0).normalized()).toRotationMatrix();
    ControlPointMatrix rigid = points*(rotation-Eigen::Matrix3d::Identity()).transpose();
    rigid.rowwise() += Eigen::RowVector3d(0.3,-0.2,0.4);
    const auto response = iga::computeUpdatedLagrangianElementResponse(space,0,rigid,material);
    require(response.internalForce.norm()<2e-13,"rigid motion must generate no internal force");
    ControlPointMatrix inverted = ControlPointMatrix::Zero(points.rows(),3);
    inverted.col(0) = -2.0*points.col(0);
    requireDomainError([&] { iga::computeUpdatedLagrangianElementResponse(space,0,inverted,material); },
                       "inverted current geometry must be rejected");
    const ControlPointMatrix collapsed = -points;
    requireDomainError([&] { iga::computeUpdatedLagrangianElementResponse(space,0,collapsed,material); },
                       "singular current geometry must be rejected");
    const ControlPointMatrix invalidReference = -points;
    requireDomainError([&] { iga::evaluateReferenceGradients(
        space.evaluateShapeFunctions(0,{0,0,0}),invalidReference,1.0); },
        "inverted reference geometry must be rejected");
}

void testGlobalAndSolverSelection()
{
    const auto space = makePatch();
    const DirectSpatialNeoHookean material;
    const auto displacement = makeDisplacements(space);
    const auto tl = iga::assembleInternalSystem(space,displacement,material);
    material.spatialCalls=0;
    const auto ul = iga::assembleInternalSystem(space,displacement,material,iga::ElementFormulation::UpdatedLagrangian);
    require(material.spatialCalls>0,"global assembly must choose UL integration");
    requireClose(ul.internalForce,tl.internalForce,3e-12,"global TL/UL internal force");
    requireClose(Eigen::MatrixXd(ul.tangentStiffness),Eigen::MatrixXd(tl.tangentStiffness),
                 3e-12,"global TL/UL stiffness");

    const auto fixed = iga::collectXiMinFaceDofs(space);
    const auto force = iga::assembleBodyForce(space,1.0,Eigen::Vector3d(0.0,0.0,-0.01));
    iga::NewtonSettings settings;
    settings.printIterations=false;
    const auto staticTL = iga::solveDeadLoad(space,material,force,fixed,settings);
    settings.elementFormulation=iga::ElementFormulation::UpdatedLagrangian;
    material.spatialCalls=0;
    const auto staticUL = iga::solveDeadLoad(space,material,force,fixed,settings);
    require(material.spatialCalls>0,"Newton settings must propagate UL selection");
    requireClose(staticUL.displacements,staticTL.displacements,2e-10,"TL/UL static solutions");

    iga::DynamicsSettings dynamic;
    dynamic.printIterations=false;
    dynamic.timeStep=0.01;
    dynamic.numberOfSteps=3;
    const Eigen::VectorXd zero=Eigen::VectorXd::Zero(iga::numDofs(space));
    const iga::TimeDependentLoad load=[force](double) { return force; };
    const auto dynamicTL=iga::solveDynamics(space,material,1.0,zero,zero,fixed,load,dynamic);
    dynamic.elementFormulation=iga::ElementFormulation::UpdatedLagrangian;
    material.spatialCalls=0;
    const auto dynamicUL=iga::solveDynamics(space,material,1.0,zero,zero,fixed,load,dynamic);
    require(material.spatialCalls>0,"dynamics settings must propagate UL selection");
    for (std::size_t i=0;i<dynamicTL.states.size();++i) {
        requireClose(dynamicUL.states[i].displacements,dynamicTL.states[i].displacements,
                     2e-11,"TL/UL dynamic displacements");
        requireClose(dynamicUL.states[i].accelerations,dynamicTL.states[i].accelerations,
                     2e-10,"TL/UL dynamic accelerations");
    }
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string,std::function<void()>>> tests{
        {"spatial tensor push-forward and shear convention",testSpatialPushForward},
        {"spatial tensor pull-back and invalid inputs",testSpatialPullBack},
        {"Neo-Hookean formulas, energy derivative and invalid inputs",testNeoHookeanMaterial},
        {"principal stretches, rotated repeated roots and tangent finite differences",testPrincipalStretchMaterial},
        {"principal-stretch invalid inputs",testPrincipalStretchInvalidInputs},
        {"UL geometry, gradients and three mappings",testGeometryAndMeasures},
        {"TL/UL equivalence and independent tangent finite difference",testElementEquivalenceAndTangent},
        {"rigid motion and invalid configurations",testRigidMotionAndInvalidStates},
        {"global assembly and static/dynamic UL selection",testGlobalAndSolverSelection}
    };
    int failures=0;
    for (const auto& test:tests) {
        try { test.second(); std::cout<<"[PASS] "<<test.first<<'\n'; }
        catch (const std::exception& error) {
            ++failures;
            std::cerr<<"[FAIL] "<<test.first<<": "<<error.what()<<'\n';
        }
    }
    if (failures) return 1;
    std::cout<<"All updated Lagrangian tests passed\n";
    return 0;
}
