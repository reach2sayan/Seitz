#pragma once

#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/math/integer_matrix.hpp>

// Operations on fractional coordinates. These mirror spglib's mat_Dmod1 /
// mat_rem1 exactly, since the unit-cell folding tolerance affects which atoms
// are judged coincident during the symmetry search.
namespace spglib::math {

// Fold a coordinate into [0, 1) with a small negative tolerance (spglib
// `mat_Dmod1`): a value within kZeroPrec below zero stays near zero rather than
// wrapping up to ~1.
[[nodiscard]] constexpr double mod1(double a) noexcept {
  double const b = a - nint(a);
  return b < -kZeroPrec ? b + 1.0 : b;
}

// Centered remainder in [-0.5, 0.5] (spglib `mat_rem1`) — the minimal-image
// offset.
[[nodiscard]] constexpr double rem1(double a) noexcept { return a - nint(a); }

[[nodiscard]] inline Vector3d mod1(Vector3d const &v) noexcept {
  return {mod1(v[0]), mod1(v[1]), mod1(v[2])};
}

[[nodiscard]] inline Vector3d rem1(Vector3d const &v) noexcept {
  return {rem1(v[0]), rem1(v[1]), rem1(v[2])};
}

// Minimal-image fractional displacement b - a, each component folded to [-0.5,
// 0.5].
[[nodiscard]] inline Vector3d frac_displacement(Vector3d const &a,
                                                Vector3d const &b) noexcept {
  return rem1(Vector3d(b - a));
}

} // namespace spglib::math
