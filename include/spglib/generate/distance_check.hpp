#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/types.hpp>

#include <optional>

namespace spglib::generate {

// Minimum-distance validity for a generated structure (PyXtal
// check_short_distances). Two atoms of types t_i and t_j may be no closer than
// `scale * (r_cov(t_i) + r_cov(t_j))`, where r_cov is the covalent radius
// (data::covalent_radius). The check covers periodic images and an atom against
// its own images, so a too-small cell is rejected as well as overlapping atoms.
struct DistanceTolerance {
  // Fraction of the summed covalent radii that sets the minimum allowed
  // separation. < 1 because a generic (non-bonded) contact in a generated
  // crystal need not reach a full bond length; the default rejects only genuine
  // clashes while leaving realistic packings valid.
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

} // namespace spglib::generate
