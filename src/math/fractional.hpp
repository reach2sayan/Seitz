#pragma once

#include <seitz/core/fractional.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

// Fractional-coordinate comparisons used by the symmetry search. The folding
// primitives themselves (wrap_to_unit_cell, nearest_offset) are public, in
// core/fractional.hpp.
namespace seitz::math {

// Minimal-image fractional displacement b - a, each component folded to
// [-0.5, 0.5].
[[nodiscard]] inline Vector3d frac_displacement(Vector3d const &a,
                                                Vector3d const &b) noexcept {
  return nearest_offset(b - a);
}

// Whether two fractional points coincide modulo lattice translations, within
// `tol` per component.
[[nodiscard]] inline bool same_point(Vector3d const &a, Vector3d const &b,
                                     double tol) noexcept {
  return approx_zero(frac_displacement(a, b), tol);
}

} // namespace seitz::math
