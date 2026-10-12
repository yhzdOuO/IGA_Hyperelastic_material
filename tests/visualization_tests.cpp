#include "visualization/BoundarySampler.h"
#include "visualization/BoundarySamplingCache.h"
#include "visualization/CantileverExample.h"

#include "mechanics/SaintVenantKirchhoffMaterial.h"

#include <Eigen/Geometry>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> openUniformKnots(int degree, int elements)
{
    std::vector<double> knots(static_cast<std::size_t>(degree + 1), 0.0);
    for (int i = 1; i < elements; ++i) knots.push_back(1.0 * i / elements);
    knots.insert(knots.end(), static_cast<std::size_t>(degree + 1), 1.0);
    return knots;
}

std::vector<double> greville(const std::vector<double>& knots, int degree)
{
    const int count = static_cast<int>(knots.size()) - degree - 1;
    std::vector<double> result(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
        for (int j = 1; j <= degree; ++j) result[i] += knots[i + j];
        result[i] /= degree;
    }
    return result;
}

nurbs::TensorProductNurbsSpace makeBox()
{
    const std::array<int, 3> degrees{2, 2, 2};
    std::array<std::vector<double>, 3> knots{
        openUniformKnots(2, 2), openUniformKnots(2, 1), openUniformKnots(2, 1)};
    std::array<std::vector<double>, 3> coordinates{
        greville(knots[0], 2), greville(knots[1], 2), greville(knots[2], 2)};
    nurbs::ControlPointMatrix points(
        coordinates[0].size() * coordinates[1].size() * coordinates[2].size(), 3);
    int global = 0;
    for (double z : coordinates[2]) {
        for (double y : coordinates[1]) {
            for (double x : coordinates[0]) points.row(global++) << 2.0 * x, y, z;
        }
    }
    return nurbs::TensorProductNurbsSpace(
        {nurbs::KnotVector(std::move(knots[0])),
         nurbs::KnotVector(std::move(knots[1])),
         nurbs::KnotVector(std::move(knots[2]))},
        degrees, std::move(points), Eigen::VectorXd::Ones(global));
}

Eigen::VectorXd affineDisplacement(const nurbs::TensorProductNurbsSpace& space,
                                   double axialStrain)
{
    Eigen::VectorXd displacement = Eigen::VectorXd::Zero(3 * space.controlPoints().rows());
    for (int point = 0; point < space.controlPoints().rows(); ++point) {
        displacement(3 * point) = axialStrain * space.controlPoints()(point, 0);
    }
    return displacement;
}

void testFieldEvaluation()
{
    const nurbs::TensorProductNurbsSpace space = makeBox();
    const double strain = 0.02;
    const iga::SaintVenantKirchhoffMaterial material(3.0, 5.0);
    const visualization::FieldEvaluator evaluator(
        space, affineDisplacement(space, strain), material,
        Eigen::Vector3d::UnitX(), -Eigen::Vector3d::UnitZ());
    const visualization::FieldPoint point = evaluator.evaluate({0.31, 0.42, 0.73});
    const double f = 1.0 + strain;
    const double greenStrain = 0.5 * (f * f - 1.0);
    const double pk2Axial = (material.lambda() + 2.0 * material.mu()) * greenStrain;
    const double expectedCauchy = f * pk2Axial;
    require(std::abs(point.cauchyStress(0, 0) - expectedCauchy) < 2.0e-12,
            "PK2-to-Cauchy push-forward is incorrect");
    require(std::abs(point.scalars[static_cast<int>(visualization::ScalarField::AxialStress)]
                     - expectedCauchy) < 2.0e-12,
            "axial stress extraction is incorrect");
    require(std::abs(point.displacement(0) - strain * point.referencePosition(0)) < 1.0e-14,
            "displacement interpolation is incorrect");
}

void testBoundarySamplingAndWinding()
{
    const nurbs::TensorProductNurbsSpace space = makeBox();
    const iga::SaintVenantKirchhoffMaterial material(3.0, 5.0);
    const visualization::FieldEvaluator evaluator(
        space, Eigen::VectorXd::Zero(3 * space.controlPoints().rows()), material,
        Eigen::Vector3d::UnitX(), -Eigen::Vector3d::UnitZ());
    visualization::SurfaceMesh mesh = visualization::BoundarySampler(3).sample(space, evaluator);
    require(mesh.vertices.size() == 90, "unexpected per-span boundary vertex count");
    require(mesh.triangleIndices.size() == 240, "unexpected triangle index count");
    require(mesh.wireIndices.size() == 160, "unexpected element-boundary line count");
    for (std::size_t triangle = 0; triangle < mesh.triangleIndices.size(); triangle += 3) {
        const visualization::SurfaceVertex& a = mesh.vertices[mesh.triangleIndices[triangle]];
        const visualization::SurfaceVertex& b = mesh.vertices[mesh.triangleIndices[triangle + 1]];
        const visualization::SurfaceVertex& c = mesh.vertices[mesh.triangleIndices[triangle + 2]];
        const Eigen::Vector3d pa(a.draw.position[0], a.draw.position[1], a.draw.position[2]);
        const Eigen::Vector3d pb(b.draw.position[0], b.draw.position[1], b.draw.position[2]);
        const Eigen::Vector3d pc(c.draw.position[0], c.draw.position[1], c.draw.position[2]);
        const Eigen::Vector3d normal(a.draw.normal[0], a.draw.normal[1], a.draw.normal[2]);
        require((pb - pa).cross(pc - pa).dot(normal) > 0.0,
                "triangle winding and analytic outward normal disagree");
    }
}

void testExaggerationDoesNotChangeStress()
{
    const nurbs::TensorProductNurbsSpace space = makeBox();
    const iga::SaintVenantKirchhoffMaterial material(3.0, 5.0);
    const visualization::FieldEvaluator evaluator(
        space, affineDisplacement(space, 0.02), material,
        Eigen::Vector3d::UnitX(), -Eigen::Vector3d::UnitZ());
    visualization::SurfaceMesh mesh = visualization::BoundarySampler(2).sample(space, evaluator);
    const float scalar = mesh.vertices.front().draw.scalar;
    const float oldX = mesh.vertices.back().draw.position[0];
    mesh.updateGeometry(8.0);
    require(mesh.vertices.front().draw.scalar == scalar,
            "deformation exaggeration must not alter stress");
    require(mesh.vertices.back().draw.position[0] != oldX,
            "deformation exaggeration must update draw coordinates");
}

void testCachedSamplingMatchesDirectEvaluation()
{
    const nurbs::TensorProductNurbsSpace space = makeBox();
    const iga::SaintVenantKirchhoffMaterial material(3.0, 5.0);
    const Eigen::VectorXd displacement = affineDisplacement(space, 0.02);
    const visualization::FieldEvaluator evaluator(
        space, displacement, material, Eigen::Vector3d::UnitX(),
        -Eigen::Vector3d::UnitZ());
    const visualization::SurfaceMesh direct =
        visualization::BoundarySampler(3).sample(space, evaluator);
    const visualization::BoundarySamplingCache cache(space, 3);
    const visualization::SurfaceMesh cached = cache.makeFrame(
        displacement, cache.evaluateScalars(displacement, material));
    require(cached.vertices.size() == direct.vertices.size(),
            "cached sampler topology mismatch");
    double maximumError = 0.0;
    for (std::size_t i = 0; i < direct.vertices.size(); ++i) {
        maximumError = std::max(maximumError,
            (cached.vertices[i].displacement - direct.vertices[i].displacement).norm());
        for (std::size_t field = 0; field < cached.vertices[i].scalars.size(); ++field) {
            maximumError = std::max(maximumError,
                std::abs(cached.vertices[i].scalars[field]
                         - direct.vertices[i].scalars[field]));
        }
    }
    require(maximumError < 3e-12,
            "precomputed sampling must match direct evaluation; error="
                + std::to_string(maximumError));
}

void testCantileverPipeline()
{
    const visualization::CantileverVisualization result =
        visualization::solveCantileverVisualization(2);
    require(result.tipDisplacement < 0.0,
            "cantilever tip must move in the gravity direction");
    require(result.newtonIterations >= 4,
            "nonlinear cantilever should require multiple Newton iterations");
    require(!result.mesh.vertices.empty() && !result.mesh.triangleIndices.empty(),
            "cantilever solve must produce a drawable boundary mesh");
    require(result.mesh.scalarMinimum < 0.0 && result.mesh.scalarMaximum > 0.0,
            "cantilever axial-stress range must include tension and compression");
}

} // namespace

int main()
{
    const std::vector<std::pair<std::string, std::function<void()>>> tests{
        {"field evaluation and Cauchy stress", testFieldEvaluation},
        {"boundary sampling and winding", testBoundarySamplingAndWinding},
        {"stress independent of exaggeration", testExaggerationDoesNotChangeStress},
        {"cached sampling matches direct evaluation",
         testCachedSamplingMatchesDirectEvaluation},
        {"cantilever solve-to-mesh pipeline", testCantileverPipeline}};
    int failures = 0;
    for (const auto& test : tests) {
        try {
            test.second();
            std::cout << "[PASS] " << test.first << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[FAIL] " << test.first << ": " << error.what() << '\n';
        }
    }
    return failures == 0 ? 0 : 1;
}
