#pragma once

#include <seitz/core/types.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>

// SVD-backed subspace primitives shared by every locus-arrangement derivation
// (point groups, rod groups). All bases are orthonormal, 3 x dim; dim == 0 is
// the origin {0}. `tol` is the singular-value cutoff — callers own their
// tolerance (the derivations legitimately differ: 1e-7 for exact point-group
// subspaces, 1e-6 for rod loci built from solved affine systems).
namespace seitz::math {

// Orthonormal basis of the column space of `m` (columns with non-negligible
// singular value).
[[nodiscard]] inline Eigen::MatrixXd column_space(Eigen::MatrixXd const &m,
                                                  double tol) {
  if (m.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU);
  Eigen::VectorXd const sv = svd.singularValues();
  auto const rank =
      std::ranges::count_if(sv, [=](double x) { return x > tol; });
  return svd.matrixU().leftCols(rank);
}

// Orthonormal basis of the null space of `m` (the domain-space directions that
// `m` annihilates).
[[nodiscard]] inline Eigen::MatrixXd null_space(Eigen::MatrixXd const &m,
                                                double tol) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullV);
  Eigen::VectorXd const sv = svd.singularValues();
  auto const rank =
      std::ranges::count_if(sv, [=](double x) { return x > tol; });
  // The last `nullity` columns of V correspond to the smallest singular values.
  return svd.matrixV().rightCols(m.cols() - rank);
}

// Orthogonal projector onto the span of an orthonormal basis.
[[nodiscard]] inline Matrix3d projector(Eigen::MatrixXd const &basis) {
  if (basis.cols() == 0) {
    return Matrix3d::Zero();
  }
  return basis * basis.transpose();
}

// Orthonormal basis of the intersection of two column spaces:
// v = a x = b y  <=>  [a  -b] [x; y] = 0, mapped back through `a`.
[[nodiscard]] inline Eigen::MatrixXd
intersect_column_spaces(Eigen::MatrixXd const &a, Eigen::MatrixXd const &b,
                        double tol) {
  if (a.cols() == 0 || b.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  Eigen::MatrixXd stacked(3, a.cols() + b.cols());
  stacked.leftCols(a.cols()) = a;
  stacked.rightCols(b.cols()) = -b;
  Eigen::MatrixXd const ker = null_space(stacked, tol); // (da+db) x k
  if (ker.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  return column_space(a * ker.topRows(a.cols()), tol);
}

} // namespace seitz::math
