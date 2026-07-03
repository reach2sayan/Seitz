#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>

namespace cppcrystal::generate {

struct DistanceTolerance {
  double scale = 0.7;
  double fallback_radius = 1.0;
};

// Smallest Cartesian distance between fractional positions `a` and `b` over all
// periodic images, searched across the 3x3x3 block of neighbouring cells (exact
// minimal image even for skewed lattices, unlike a single nearest-offset fold).
// When `include_origin` is false the zero offset is skipped — use that for an
// atom against its own images (a == b), so the result is the nearest non-trivial
// self-image rather than 0. With `aperiodic_axis` set (layer cells), that axis
// is not periodic: its component is kept as the raw difference and no neighbour
// cells are searched along it.
[[nodiscard]] double
minimum_image_distance(Vector3d const &a, Vector3d const &b,
                       Matrix3d const &lattice, bool include_origin = true,
                       std::optional<int> aperiodic_axis = std::nullopt) noexcept;

// True iff every pair of atoms in `cell` (including each atom against its own
// periodic images) is at least its type-pair minimum distance apart.
[[nodiscard]] bool distances_valid(Cell const &cell,
                                   DistanceTolerance tol = {}) noexcept;

// True iff every pair of atoms in a non-periodic cluster is at least its
// type-pair minimum distance apart. `coordinates` are Cartesian (row i = atom
// i), `types[i]` the type of atom i. The 0D analogue of distances_valid: no
// lattice and no periodic images, so it is a plain all-pairs Euclidean check —
// used by the point-group cluster generator (generate::random_cluster).
[[nodiscard]] bool cluster_distances_valid(Positions const &coordinates,
                                           Types const &types,
                                           DistanceTolerance tol = {}) noexcept;

} // namespace cppcrystal::generate
