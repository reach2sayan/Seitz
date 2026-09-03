#pragma once

#include <cppcrystal/core/types.hpp>

#include <Eigen/Dense>

#include <cmath>
#include <optional>

// Small numerical helpers backing the symmetry search, where exact semantics
// matter: integer rounding, a checked real inverse, and an exact inverse for
// the unimodular integer rotation matrices. General matrix arithmetic uses
// Eigen operators at the call site.
namespace cppcrystal::math {

// Round to the nearest integer, ties away from zero (matching the truncation
// behaviour the symmetry search depends on).
[[nodiscard]] constexpr int nint(double a) noexcept {
  return a < 0.0 ? static_cast<int>(a - 0.5) : static_cast<int>(a + 0.5);
}

[[nodiscard]] inline auto round_to_int(MatrixExpr auto const &a) noexcept {
  return a.unaryExpr([](double x) { return nint(x); }).eval();
}

// Displacement from each coordinate to its nearest integer, x - nint(x). For a
// fractional position this is the offset to the nearest lattice point, the
// quantity that is invariant under integer lattice translations.
[[nodiscard]] inline auto nearest_offset(MatrixExpr auto const &a) noexcept {
  return a.unaryExpr([](double x) { return x - nint(x); });
}

// True iff every entry of `a` is within `symprec` of an integer.
[[nodiscard]] inline bool is_int_matrix(const Matrix3d &a,
                                        double symprec) noexcept {
  return a.unaryExpr([](double x) {
            return std::abs(nint(x) - x);
          }).maxCoeff() <= symprec;
}

// Checked real inverse: std::nullopt when |det| < precision.
[[nodiscard]] inline std::optional<Matrix3d>
inverse(Matrix3d const &a, double precision) noexcept {
  double const det = a.determinant();
  if (std::abs(det) < precision) {
    return std::nullopt;
  }
  return a.inverse();
}

// Exact integer inverse of a unimodular matrix (|det| == 1); std::nullopt
// otherwise. Lets us invert integer rotation matrices with no floating-point
// round-trip.
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

} // namespace cppcrystal::math
