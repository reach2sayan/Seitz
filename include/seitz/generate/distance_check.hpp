#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>

#pragma GCC visibility push(default)

namespace seitz::generate {

struct DistanceTolerance {
  double scale = 0.7;
  double fallback_radius = 1.0;
};

// Smallest Cartesian distance between fractional positions `a` and `b` over
// (a) periodic images: each periodic component is folded to [-0.5, 0.5] and
// the neighbouring cells are searched along the periodic axes
// (b) aperiodic images : keeps raw difference with no neighbour search.
//
// `Images::nontrivial` skips the zero offset —
//  use it for an atom against its own images (a == b),
enum class Images { all, nontrivial };

[[nodiscard]] double minimum_image_distance(
    Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
    CellPeriodicity const &periodicity, Images images = Images::all) noexcept;

// True iff every pair of atoms in `cell` (including each atom against its own
// periodic images) is at least its type-pair minimum distance apart, under the
// cell's own periodicity.
[[nodiscard]] bool distances_valid(Cell const &cell,
                                   DistanceTolerance tol = {});

} // namespace seitz::generate

#pragma GCC visibility pop
