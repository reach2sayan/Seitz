#pragma once

#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

// Operations on fractional coordinates. The unit-cell folding tolerance affects
// which atoms are judged coincident during the symmetry search. The
// minimal-image offset x - nint(x) lives in integer_matrix.hpp as
// `nearest_offset`.
namespace cppcrystal::math {

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

// Minimal-image fractional displacement b - a, each component folded to [-0.5,
// 0.5].
[[nodiscard]] inline Vector3d frac_displacement(Vector3d const &a,
                                                Vector3d const &b) noexcept {
  return nearest_offset(b - a);
}

// Whether two fractional points coincide modulo lattice translations, within
// `tol` per component.
[[nodiscard]] inline bool same_point(Vector3d const &a, Vector3d const &b,
                                     double tol) noexcept {
  return frac_displacement(a, b).cwiseAbs().maxCoeff() < tol;
}

} // namespace cppcrystal::math
