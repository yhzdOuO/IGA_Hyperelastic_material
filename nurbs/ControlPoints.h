#ifndef CONTROLPOINTS_H
#define CONTROLPOINTS_H

#pragma once
#include <Eigen/Core>
#include <stdexcept>

namespace nurbs {

// -----------------------------------------------------------------------
// Three-dimensional control-point storage used by the solid solver.  Geometry
// and weights deliberately live in separate contiguous Eigen objects so that
// Jacobians can be formed directly as P.transpose() * dN_dxi.
// -----------------------------------------------------------------------
using ControlPointMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3>;

struct ControlPointData {
    ControlPointMatrix P;
    Eigen::VectorXd weights;

    /** @brief Returns the number of control points. @return Control-point count. */
    int count() const { return static_cast<int>(P.rows()); }

    /** @brief Returns the highest valid control-point index. @return Last control-point index. */
    int lastIndex() const { return count() - 1; }

    /** @brief Validates matching sizes, finite coordinates, and positive finite weights. */
    void validate() const
    {
        if (P.rows() != weights.size()) {
            throw std::invalid_argument("ControlPointData: coordinate and weight counts differ");
        }
        if (!P.allFinite() || !weights.allFinite() || (weights.array() <= 0.0).any()) {
            throw std::invalid_argument("ControlPointData: coordinates must be finite and weights positive");
        }
    }
};

/**
 * @brief Checks the B-spline relation m = n+p+1.
 * @param numKnots Number of knots, m+1.
 * @param p B-spline degree.
 * @param numCtrlPts Number of control points, n+1.
 * @return True when the knot and control-point counts are consistent.
 */
inline bool isKnotControlConsistent(int numKnots, int p, int numCtrlPts)
{
    const int m = numKnots - 1;
    const int n = numCtrlPts - 1;
    return m == n + p + 1;
}

} // namespace nurbs

#endif // CONTROLPOINTS_H
