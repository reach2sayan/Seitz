#pragma once

#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

// The handful of numeric primitives the public headers need inline: integer
// rounding and the two fractional-coordinate foldings. The rest of the
// numerics (checked/exact inverses, integer-matrix predicates, subspaces) is
// private to src/math.
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

// Wrap a coordinate back into the unit cell [0, 1) with a small negative
// tolerance: a value within kZeroPrec below zero stays near zero rather than
// wrapping up to ~1.
[[nodiscard]] constexpr double wrap_to_unit_cell(double a) noexcept {
  double const b = a - nint(a);
  return b < -kZeroPrec ? b + 1.0 : b;
}

[[nodiscard]] inline Vector3d wrap_to_unit_cell(Vector3d const &v) noexcept {
  return {wrap_to_unit_cell(v[0]), wrap_to_unit_cell(v[1]),
          wrap_to_unit_cell(v[2])};
}

} // namespace cppcrystal::math
