#ifndef TENSORPRODUCTNURBSSPACE_H
#define TENSORPRODUCTNURBSSPACE_H

#include "ControlPoints.h"
#include "KnotVector.h"

#include <Eigen/Core>
#include <array>
#include <vector>

namespace nurbs {

struct ShapeFunctionData {
    /// Rational NURBS basis values R_A at one integration point.
    Eigen::VectorXd values;
    /// Rows contain [dR_A/dxi, dR_A/deta, dR_A/dzeta].
    Eigen::Matrix<double, Eigen::Dynamic, 3> parametricGradients;
};

/** @brief Active basis data and connectivity at an arbitrary parameter point. */
struct PointShapeFunctionData {
    std::vector<int> connectivity; // &connectivity[local_index] = global_index
    ShapeFunctionData shape;
};

/** @brief Stores a three-dimensional tensor-product NURBS patch and its L0 tables. */
class TensorProductNurbsSpace {
public:
    /**
     * @brief Builds and validates a NURBS solid space, span tables, and Gauss basis tables.
     * @param knotVectors Knot vectors in the xi, eta, and zeta directions.
     * @param degrees Polynomial degrees in the three parameter directions.
     * @param controlPoints Global control-point coordinates; xi index varies fastest.
     * @param weights Positive global NURBS weights in the same order as controlPoints.
     */
    TensorProductNurbsSpace(std::array<KnotVector, 3> knotVectors,
                           std::array<int, 3> degrees,
                           ControlPointMatrix controlPoints,
                           Eigen::VectorXd weights);

    /** @brief Returns a directional knot vector. @param direction 0=xi, 1=eta, 2=zeta. @return Requested knot vector. */
    const KnotVector& knotVector(int direction) const;

    /** @brief Returns a directional polynomial degree. @param direction 0=xi, 1=eta, 2=zeta. @return Polynomial degree. */
    int degree(int direction) const;

    /** @brief Returns the directional control-point count. @param direction Parameter direction. @return Control-point count. */
    int numControlPoints(int direction) const;

    /** @brief Returns the number of nonzero spans in one direction. @param direction Parameter direction. @return Directional element count. */
    int numElements(int direction) const;

    /** @brief Returns the total tensor-product element count. @return Total element count. */
    int numElements() const;

    /** @brief Returns the number of active basis functions per element. @return Local basis count n_en. */
    int numElementNodes() const;

    /** @brief Returns the Gauss-point count in one direction. @param direction Parameter direction. @return Directional Gauss-point count. */
    int numGaussPoints(int direction) const;

    /**
     * @brief Returns one Gauss-Legendre abscissa on the parent interval [-1,1].
     * @param direction Parameter direction whose degree determines the rule order.
     * @param gauss Directional Gauss-point index.
     * @return Parent-coordinate Gauss point before mapping to a knot span.
     */
    double gaussPoint(int direction, int gauss) const;

    /**
     * @brief Returns one directional Gauss-Legendre weight on [-1,1].
     * @param direction Parameter direction whose degree determines the rule order.
     * @param gauss Directional Gauss-point index.
     * @return One-dimensional quadrature weight.
     */
    double gaussWeight(int direction, int gauss) const;

    /**
     * @brief Returns the knot-span index of a directional element.
     * @param direction Parameter direction.
     * @param directionalElement Zero-based element index in that direction.
     * @return Knot-span index compatible with FindSpan().
     */
    int span(int direction, int directionalElement) const;

    /**
     * @brief Unflattens a global element number with xi varying fastest.
     * @param element Flattened element index.
     * @return Directional element indices {e_xi, e_eta, e_zeta}.
     */
    std::array<int, 3> elementIndices(int element) const;

    /** @brief Returns the three knot-span indices for an element. @param element Element index. @return {s_xi, s_eta, s_zeta}. */
    std::array<int, 3> elementSpans(int element) const;

    /**
     * @brief Evaluates the lazy local-to-global basis mapping.
     * @param localBasis Local basis index with xi varying fastest.
     * @param element Element index.
     * @return Global control-point/basis index.
     */
    int IEN(int localBasis, int element) const;

    /** @brief Generates all local-to-global basis indices for one element. @param element Element index. @return Connectivity in local basis order. */
    std::vector<int> elementConnectivity(int element) const;

    /** @brief Gathers an element's control-point coordinates. @param element Element index. @return n_en by 3 coordinate matrix. */
    ControlPointMatrix elementControlPoints(int element) const;

    /** @brief Gathers an element's NURBS weights. @param element Element index. @return Local weight vector. */
    Eigen::VectorXd elementWeights(int element) const;

    /**
     * @brief Evaluates rational shape values and parameter gradients from precomputed tables.
     * @param element Element index.
     * @param gaussIndices Directional Gauss-point indices {q_xi, q_eta, q_zeta}.
     * @return Shape values and their xi/eta/zeta derivatives.
     */
    ShapeFunctionData evaluateShapeFunctions(
        int element, const std::array<int, 3>& gaussIndices) const;

    /**
     * @brief Evaluates active NURBS functions at an arbitrary parameter point.
     * @param parameters Coordinates {xi, eta, zeta} inside the patch domain.
     * @return Active connectivity, shape values, and parameter gradients.
     */
    PointShapeFunctionData evaluateShapeFunctionsAt(
        const std::array<double, 3>& parameters) const;

    /** @brief Returns det(J_parent) for the parent-to-parameter map. @param element Element index. @return Parent Jacobian determinant. */
    double parentJacobianDeterminant(int element) const;

    /** @brief Returns the three-dimensional tensor-product Gauss weight. @param gaussIndices Directional Gauss indices. @return Product of the three univariate weights. */
    double gaussWeight(const std::array<int, 3>& gaussIndices) const;

    /** @brief Returns the global control-point coordinate matrix. @return Read-only coordinate matrix reference. */
    const ControlPointMatrix& controlPoints() const { return controlPoints_; }

    /** @brief Returns the global NURBS weight vector. @return Read-only weight-vector reference. */
    const Eigen::VectorXd& weights() const { return weights_; }

private:
    struct GaussRule {
        std::vector<double> points;
        std::vector<double> weights;
    };

    struct UnivariateTable {
        int numElements = 0;
        int numGauss = 0;
        int numBasis = 0;
        std::vector<double> data;

        /**
         * @brief Reads one entry from the flattened univariate basis table.
         * @param element Directional element index.
         * @param gauss Directional Gauss-point index.
         * @param derivative Derivative order stored in the table (0 or 1).
         * @param basis Local univariate basis index.
         * @return Requested basis value or first derivative.
         */
        double value(int element, int gauss, int derivative, int basis) const;
    };

    /** @brief Builds an order-point Gauss-Legendre rule on [-1,1]. @param order Rule order. @return Gauss points and weights. */
    static GaussRule makeGaussLegendreRule(int order);

    /** @brief Scans the knot vectors and stores all nonzero span indices. */
    void buildSpans();

    /** @brief Precomputes univariate basis values and first derivatives at Gauss points. */
    void buildTables();

    /** @brief Validates degrees, control-net dimensions, coordinates, and weights. */
    void validate() const;

    std::array<KnotVector, 3> knotVectors_;
    std::array<int, 3> degrees_;
    std::array<int, 3> controlPointCounts_{};
    std::array<std::vector<int>, 3> spans_;
    std::array<GaussRule, 3> gaussRules_;
    std::array<UnivariateTable, 3> tables_;
    ControlPointMatrix controlPoints_;
    Eigen::VectorXd weights_;
};

} // namespace nurbs

#endif // TENSORPRODUCTNURBSSPACE_H
