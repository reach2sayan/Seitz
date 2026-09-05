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

// Smallest Cartesian distance between fractional `a` and `b`: periodic
// components folded to [-0.5, 0.5] with a neighbour search along those axes,
// aperiodic ones kept as the raw difference. `Images::nontrivial` skips the
// zero offset -- for an atom against its own images (a == b).
enum class Images { all, nontrivial };

[[nodiscard]] double minimum_image_distance(
    Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
    CellPeriodicity const &periodicity, Images images = Images::all) noexcept;

// True iff every atom pair in `cell` (each atom against its own images
// included) clears its type-pair minimum distance, under the cell's
// periodicity.
[[nodiscard]] bool distances_valid(Cell const &cell,
                                   DistanceTolerance tol = {});

} // namespace seitz::generate

#pragma GCC visibility pop
