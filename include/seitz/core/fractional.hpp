#pragma once

#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <cmath>

// The handful of numeric primitives the public headers need inline: integer
// rounding and the two fractional-coordinate foldings. The rest of the
// numerics (checked/exact inverses, integer-matrix predicates, subspaces) is
// private to src/math.

#pragma GCC visibility push(default)

namespace seitz::math {

// Round to the nearest integer, ties away from zero (matching the truncation
// behaviour the symmetry search depends on).
//
// Branchless on purpose. This is the innermost operation in the library --
// nearest_offset, wrap_to_unit_cell, minimal_image and coincident all apply it
// elementwise -- and the ternary form it replaces compiled to a compare and a
// jump that blocked vectorization; the copysign form compiles to andpd/orpd
// and lets GCC do four lanes at a time (cvttpd2dq). Bit-identical to the
// ternary for every double, -0.0 included: copysign(0.5, -0.0) is -0.5, and
// the conversion truncates toward zero either way.
//
// std::copysign is only constexpr from C++23 (P0533) and MSVC has not
// implemented that yet, so constant evaluation takes the ternary form -- which
// agrees with copysign on every double for this use, including -0.0, where
// both -0.5 and +0.5 truncate to 0. Runtime keeps the vectorizable form.
[[nodiscard]] constexpr int nint(double a) noexcept {
  if consteval {
    return static_cast<int>(a + (a < 0.0 ? -0.5 : 0.5));
  }
  return static_cast<int>(a + std::copysign(0.5, a));
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

} // namespace seitz::math

#pragma GCC visibility pop
