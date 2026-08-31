#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

namespace cppcrystal::generate {

struct DistanceTolerance {
  double scale = 0.7;
  double fallback_radius = 1.0;
};

// Smallest Cartesian distance between fractional positions `a` and `b` over
// all periodic images: each periodic component is folded to [-0.5, 0.5] and
// the neighbouring cells are searched along the periodic axes (exact minimal
// image even for skewed lattices); each aperiodic component keeps its raw
// difference with no neighbour search. This covers every family with one
// function — 3D (all periodic), layer (one aperiodic), rod (one periodic) and
// cluster (none periodic — a plain Euclidean distance). When `include_origin`
// is false the zero offset is skipped — use that for an atom against its own
// images (a == b), so the result is the nearest non-trivial self-image rather
// than 0 (infinite when no axis is periodic).
[[nodiscard]] double
minimum_image_distance(Vector3d const &a, Vector3d const &b,
                       Matrix3d const &lattice,
                       CellPeriodicity const &periodicity,
                       bool include_origin = true) noexcept;

// True iff every pair of atoms in `cell` (including each atom against its own
// periodic images) is at least its type-pair minimum distance apart, under
// the cell's own periodicity.
[[nodiscard]] bool distances_valid(Cell const &cell,
                                   DistanceTolerance tol = {}) noexcept;

// True iff every pair of atoms in a non-periodic cluster is at least its
// type-pair minimum distance apart. `coordinates` are Cartesian (row i = atom
// i), `types[i]` the type of atom i — used by generate::random_cluster.
[[nodiscard]] bool cluster_distances_valid(Positions const &coordinates,
                                           Types const &types,
                                           DistanceTolerance tol = {}) noexcept;

} // namespace cppcrystal::generate
