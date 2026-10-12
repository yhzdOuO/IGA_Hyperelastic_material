#ifndef KNOTVECTOR_H
#define KNOTVECTOR_H

#pragma once
#include <vector>
#include <cstddef>
#include <stdexcept>
#include <algorithm>

namespace nurbs {

// -----------------------------------------------------------------------
// KnotVector
//
// Represents U = { U_0, U_1, ..., U_m }, a nondecreasing sequence of real
// numbers, as defined in Piegl & Tiller, "The NURBS Book" (2nd ed.), Ch. 2.
//
// Relationship with degree p and number of control points (n+1):
//     m = n + p + 1        (m+1 knots, n+1 control points, degree p)
// -----------------------------------------------------------------------
class KnotVector {
public:
    /** @brief Creates an empty knot vector for deferred initialization. */
    KnotVector() = default;

    /**
     * @brief Creates a validated, nondecreasing knot vector.
     * @param knots Knot values in nondecreasing order; at least two values are required.
     */
    explicit KnotVector(std::vector<double> knots)
        : knots_(std::move(knots))
    {
        if (knots_.size() < 2) {
            throw std::invalid_argument("KnotVector: need at least 2 knots");
        }
        if (!std::is_sorted(knots_.begin(), knots_.end())) {
            throw std::invalid_argument("KnotVector: knots must be nondecreasing");
        }
    }

    /** @brief Returns knot U[i]. @param i Zero-based knot index. @return Knot value. */
    double  operator[](int i) const { return knots_[static_cast<std::size_t>(i)]; }

    /** @brief Returns mutable knot U[i]. @param i Zero-based knot index. @return Knot reference. */
    double& operator[](int i)       { return knots_[static_cast<std::size_t>(i)]; }

    /** @brief Returns the number of knots, m+1. @return Knot count. */
    int size() const { return static_cast<int>(knots_.size()); }

    /** @brief Returns the highest valid knot index, m. @return Last knot index. */
    int lastIndex() const { return size() - 1; }

    /**
     * @brief Infers the highest control-point index from m = n+p+1.
     * @param p B-spline degree.
     * @return Highest control-point index n.
     */
    int inferLastCtrlPtIndex(int p) const { return lastIndex() - p - 1; }

    /** @brief Returns the first knot value. @return First knot. */
    double frontValue() const { return knots_.front(); }

    /** @brief Returns the last knot value. @return Last knot. */
    double backValue()  const { return knots_.back(); }

    /** @brief Returns read-only access to the underlying contiguous knot storage. @return Knot storage reference. */
    const std::vector<double>& raw() const { return knots_; }

private:
    std::vector<double> knots_;
};

} // namespace nurbs

#endif // KNOTVECTOR_H
