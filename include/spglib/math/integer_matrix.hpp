#pragma once

#include <spglib/core/types.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <optional>

// Small numerical helpers that survive the move from C to Eigen. The bulk of
// spglib's mathfunc.c (matrix multiply/add/transpose/determinant/trace) is
// gone, replaced by Eigen operators at the call site. What remains here are the
// few routines whose *exact* semantics the symmetry search depends on: spglib's
// integer rounding, its checked real inverse, and an exact inverse for the
// unimodular integer rotation matrices.
namespace spglib::math {

// Round to the nearest integer, ties away from zero. Exact port of spglib's
// `mat_Nint`, including its truncation behaviour, so detection matches the
// reference.
[[nodiscard]] constexpr int nint(double a) noexcept {
  return a < 0.0 ? static_cast<int>(a - 0.5) : static_cast<int>(a + 0.5);
}

template <typename Derived>
[[nodiscard]] inline auto
round_to_int(Eigen::MatrixBase<Derived> const &a) noexcept {
  return a.unaryExpr([](double x) { return nint(x); }).eval();
}

// True iff every entry of `a` is within `symprec` of an integer
// (spglib `mat_is_int_matrix`).
[[nodiscard]] inline bool is_int_matrix(const Matrix3d &a,
                                        double symprec) noexcept {
  return a.unaryExpr([](double x) {
            return std::abs(nint(x) - x);
          }).maxCoeff() <= symprec;
}

// Metric (Gram) tensor L^T . L for a lattice whose columns are the basis
// vectors (spglib `mat_get_metric`).
[[nodiscard]] inline Matrix3d metric_tensor(Matrix3d const &lattice) noexcept {
  return lattice.transpose() * lattice;
}

// Checked real inverse: std::nullopt when |det| < precision
// (spglib `mat_inverse_matrix_d3`).
[[nodiscard]] inline std::optional<Matrix3d>
inverse(Matrix3d const &a, double precision) noexcept {
  double const det = a.determinant();
  if (std::abs(det) < precision) {
    return std::nullopt;
  }
  return a.inverse();
}

// b^-1 . a . b (spglib `mat_get_similar_matrix_d3`); nullopt when b is
// singular.
[[nodiscard]] inline std::optional<Matrix3d>
similar(Matrix3d const &a, Matrix3d const &b, double precision) noexcept {
  auto const binv = inverse(b, precision);
  if (!binv) {
    return std::nullopt;
  }
  return *binv * a * b;
}

// Exact integer inverse of a unimodular matrix (|det| == 1); std::nullopt
// otherwise. Lets us invert integer rotation matrices with no floating-point
// round-trip. The adjugate entries match spglib's mat_inverse_matrix_d3 layout.
[[nodiscard]] inline std::optional<Matrix3i>
integer_inverse(Matrix3i const &a) noexcept {
  int const det = a.determinant();
  if (det != 1 && det != -1) {
    return std::nullopt;
  }
  Matrix3i adj; // adjugate = transpose of the cofactor matrix
  adj(0, 0) = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1);
  adj(0, 1) = a(0, 2) * a(2, 1) - a(0, 1) * a(2, 2);
  adj(0, 2) = a(0, 1) * a(1, 2) - a(0, 2) * a(1, 1);
  adj(1, 0) = a(1, 2) * a(2, 0) - a(1, 0) * a(2, 2);
  adj(1, 1) = a(0, 0) * a(2, 2) - a(0, 2) * a(2, 0);
  adj(1, 2) = a(0, 2) * a(1, 0) - a(0, 0) * a(1, 2);
  adj(2, 0) = a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0);
  adj(2, 1) = a(0, 1) * a(2, 0) - a(0, 0) * a(2, 1);
  adj(2, 2) = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
  return det == 1 ? adj : Matrix3i(-adj);
}

} // namespace spglib::math
