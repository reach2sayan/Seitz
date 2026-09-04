#pragma once

#include <cppcrystal/core/fractional.hpp>
#include <cppcrystal/core/types.hpp>

#include <Eigen/Dense>

#include <array>
#include <cmath>
#include <cstddef>
#include <optional>

// Numerical helpers backing the symmetry search where exact semantics matter:
// a checked real inverse, and the exact integer adjugate/determinant that give
// the unimodular rotation matrices an inverse with no floating-point round
// trip. The rounding primitives live in core/fractional.hpp, which the public
// headers need inline. General matrix arithmetic uses Eigen operators at the
// call site.
namespace cppcrystal::math {

// Row-major 3x3 integer proxy. Eigen 5's fixed-size matrices are literal types
// for construction and access, but not for arithmetic, so anything that has to
// run at compile time works on this instead — see core/centering.hpp, which
// derives its inverse tables through adjugate() at consteval time.
using Mat3Rows = std::array<int, 9>;
using Row3 = std::array<int, 3>;

[[nodiscard]] constexpr Row3 row_of(Mat3Rows const &m, std::size_t i) noexcept {
  return {m[i * 3], m[i * 3 + 1], m[i * 3 + 2]};
}

[[nodiscard]] constexpr Row3 cross(Row3 const &u, Row3 const &v) noexcept {
  return {u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
          u[0] * v[1] - u[1] * v[0]};
}

// Scalar triple product row0 . (row1 x row2).
[[nodiscard]] constexpr int determinant(Mat3Rows const &m) noexcept {
  Row3 const r0 = row_of(m, 0);
  Row3 const c = cross(row_of(m, 1), row_of(m, 2));
  return r0[0] * c[0] + r0[1] * c[1] + r0[2] * c[2];
}

// Adjugate: the transpose of the cofactor matrix, so the cofactor rows (each a
// cross product of two rows of m) land in the COLUMNS of the result. Exact in
// integer arithmetic, and M . adjugate(M) == determinant(M) . I identically.
[[nodiscard]] constexpr Mat3Rows adjugate(Mat3Rows const &m) noexcept {
  Row3 const c0 = cross(row_of(m, 1), row_of(m, 2));
  Row3 const c1 = cross(row_of(m, 2), row_of(m, 0));
  Row3 const c2 = cross(row_of(m, 0), row_of(m, 1));
  return {c0[0], c1[0], c2[0], //
          c0[1], c1[1], c2[1], //
          c0[2], c1[2], c2[2]};
}

// Row-major view of an Mat3Rows as an Eigen matrix, and its inverse mapping.
using RowMajor3i = Eigen::Matrix<int, 3, 3, Eigen::RowMajor>;

[[nodiscard]] inline Eigen::Map<RowMajor3i const>
as_matrix(Mat3Rows const &m) noexcept {
  return Eigen::Map<RowMajor3i const>(m.data());
}

[[nodiscard]] inline Mat3Rows as_rows(Matrix3i const &m) noexcept {
  Mat3Rows rows{};
  Eigen::Map<RowMajor3i>(rows.data()) = m;
  return rows;
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
// otherwise. For |det| == 1 the inverse *is* the adjugate up to sign, so this
// stays in integer arithmetic with no floating-point round trip.
[[nodiscard]] inline std::optional<Matrix3i>
integer_inverse(Matrix3i const &a) noexcept {
  Mat3Rows const rows = as_rows(a);
  int const det = determinant(rows);
  if (det != 1 && det != -1) {
    return std::nullopt;
  }
  Matrix3i const adj = as_matrix(adjugate(rows));
  return det == 1 ? adj : Matrix3i(-adj);
}

} // namespace cppcrystal::math
