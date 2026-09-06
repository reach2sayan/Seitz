#pragma once

#include <seitz/core/fractional.hpp>
#include <seitz/core/types.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <ranges>
#include <vector>

// Where the symmetry search needs exact semantics: a checked real inverse, and
// the integer adjugate/determinant giving a unimodular rotation its inverse
// with no floating-point round trip. Rounding lives in core/fractional.hpp,
// which the public headers need inline; general arithmetic is Eigen at the call
// site.
namespace seitz::math {

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

// The |det P| coset representatives of Z^3 / P Z^3, as fractional points of the
// P-cell: P^-1 n folded into [0, 1)^3 for every integer n of the frame
// surrounding P [0, 1)^3. Kept in integer arithmetic -- P^-1 n = adjugate(P) n
// / det P, reduced mod det -- so coinciding cosets compare exactly and no
// tolerance is involved. Empty for a singular P.
[[nodiscard]] inline std::vector<Vector3d>
lattice_points_in_cell(Matrix3i const &p) {
  Mat3Rows const rows = as_rows(p);
  int const det = determinant(rows);
  if (det == 0) {
    return {};
  }
  int const modulus = std::abs(det);
  Matrix3i const adj = as_matrix(adjugate(rows));
  // Bounding box of the parallelepiped: each corner is an independent 0/1
  // choice per column, so a row's extent is its positive minus its negative
  // entry sum.
  Vector3i const frame =
      p.cwiseMax(0).rowwise().sum() - p.cwiseMin(0).rowwise().sum();

  auto const residue = [&](int numerator) {
    int const r = (det > 0 ? numerator : -numerator) % modulus;
    return r < 0 ? r + modulus : r;
  };
  std::vector<std::array<int, 3>> keys{
      std::from_range,
      std::views::cartesian_product(std::views::iota(0, frame[0]),
                                    std::views::iota(0, frame[1]),
                                    std::views::iota(0, frame[2])) |
          std::views::transform([&](auto const &n) {
            auto const [i, j, k] = n;
            Vector3i const m = adj * Vector3i(i, j, k);
            return std::array{residue(m[0]), residue(m[1]), residue(m[2])};
          })};
  std::ranges::sort(keys);
  auto const [first, last] = std::ranges::unique(keys);
  keys.erase(first, last);
  return {std::from_range,
          keys | std::views::transform([&](std::array<int, 3> const &key) {
            return Vector3d(Vector3d(key[0], key[1], key[2]) / modulus);
          })};
}

} // namespace seitz::math
