#pragma once

#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <cmath>

// The numeric primitives the public headers need inline: integer rounding and
// the two fractional foldings. The rest (checked/exact inverses, integer-matrix
// predicates, subspaces) is private to src/math.

#pragma GCC visibility push(default)

namespace seitz::math {

// Round to nearest, ties away from zero. Branchless on purpose: this is the
// library's innermost operation (nearest_offset, wrap_to_unit_cell,
// minimal_image, coincident all apply it elementwise) and the ternary form
// blocked vectorization, where copysign gives andpd/orpd and four lanes of
// cvttpd2dq. Bit-identical to the ternary for every double, -0.0 included.
[[nodiscard]] constexpr int nint(double a) noexcept {
  if consteval {
    return static_cast<int>(a + (a < 0.0 ? -0.5 : 0.5));
  }
  return static_cast<int>(a + std::copysign(0.5, a));
}

[[nodiscard]] inline auto round_to_int(MatrixExpr auto const &a) noexcept {
  return a.unaryExpr([](double x) { return nint(x); }).eval();
}

// x - nint(x) per coordinate: the offset to the nearest lattice point, the
// quantity invariant under integer lattice translations.
[[nodiscard]] inline auto nearest_offset(MatrixExpr auto const &a) noexcept {
  return a.unaryExpr([](double x) { return x - nint(x); });
}

// Into [0, 1), with a negative tolerance: a value within kZeroPrec below zero
// stays near zero instead of wrapping to ~1.
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
