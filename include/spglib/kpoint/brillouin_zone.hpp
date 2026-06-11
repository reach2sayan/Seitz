#pragma once

#include <spglib/core/types.hpp>

#include <cstddef>
#include <optional>
#include <vector>

// Brillouin-zone relocation of a reciprocal sampling grid. Port of the
// BZ-related functions of spglib's kpoint.c. Each grid point is moved to the
// reciprocal-lattice image closest to the origin; points on a BZ boundary
// (multiple equidistant images within tolerance) are duplicated. The map is
// indexed on the doubled mesh (bzmesh = 2*mesh).
namespace spglib::kpoint {

// The relocated BZ grid: `bz_grid_address` holds the in-BZ grid points (the
// first prod(mesh) are the closest image of each input point at its original
// index, followed by the boundary duplicates); `bz_map[bzgp]` maps a
// doubled-mesh grid-point index to its `bz_grid_address` index, or std::nullopt
// where no BZ point falls (replacing spglib's `num_bzmesh` sentinel).
struct BzGrid {
  std::vector<Vector3i> bz_grid_address;
  std::vector<std::optional<std::size_t>> bz_map;
};

// Relocate `grid_address` (the reducible grid of `mesh`) into the first
// Brillouin zone defined by `rec_lattice` (columns = reciprocal basis vectors).
// Port of kpt_relocate_dense_BZ_grid_address.
[[nodiscard]] BzGrid
relocate_BZ_grid_address(std::vector<Vector3i> const &grid_address,
                         Vector3i const &mesh, Matrix3d const &rec_lattice,
                         Vector3i const &is_shift);

// The BZ grid-point indices that one grid address maps to under each reciprocal
// rotation: rotate on the doubled mesh, then look the result up in `bz_map`
// (std::nullopt where unmapped). Port of
// kpt_get_dense_BZ_grid_points_by_rotations.
[[nodiscard]] std::vector<std::optional<std::size_t>>
BZ_grid_points_by_rotations(
    Vector3i const &address_orig, std::vector<Matrix3i> const &rot_reciprocal,
    Vector3i const &mesh, Vector3i const &is_shift,
    std::vector<std::optional<std::size_t>> const &bz_map);

} // namespace spglib::kpoint
