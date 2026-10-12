#ifndef NURBSALGORITHM_H
#define NURBSALGORITHM_H

#include "KnotVector.h"

#include <Eigen/Core>
#include <vector>

namespace nurbs {

/**
 * @brief Finds the knot span containing a parameter using Algorithm A2.1.
 * @param n Index of the last control point; the control-point count is n+1.
 * @param p B-spline degree.
 * @param u Parameter value in the valid interval [U[p], U[n+1]].
 * @param U Knot vector associated with the basis.
 * @return Span index i such that u is in [U[i], U[i+1]); returns n at the right endpoint.
 */
int FindSpan(int n, int p, double u, const KnotVector& U);

/**
 * @brief Evaluates the p+1 nonzero B-spline basis functions using Algorithm A2.2.
 * @param i Knot span index, normally obtained from FindSpan().
 * @param u Parameter value at which the basis is evaluated.
 * @param p B-spline degree.
 * @param U Knot vector associated with the basis.
 * @param N Output values, where N[j] = N_(i-p+j,p)(u).
 */
void BasisFuns(int i, double u, int p, const KnotVector& U,
               std::vector<double>& N);

/**
 * @brief Evaluates nonzero B-spline basis functions and derivatives using Algorithm A2.3.
 * @param i Knot span index, normally obtained from FindSpan().
 * @param u Parameter value at which the basis is evaluated.
 * @param p B-spline degree.
 * @param nDer Highest derivative order requested; orders above p are returned as zero.
 * @param U Knot vector associated with the basis.
 * @param ders Output matrix with ders(k,j) = d^k N_(i-p+j,p)/du^k.
 */
void DersBasisFuns(int i, double u, int p, int nDer, const KnotVector& U,
                   Eigen::MatrixXd& ders);

} // namespace nurbs

#endif // NURBSALGORITHM_H
