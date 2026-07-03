#pragma once

#include <cppcrystal/core/types.hpp>

#include <cstddef>
#include <vector>

// Reciprocal-space sampling grid: the mapping between integer grid addresses and
// linear grid-point indices (3D path). `is_shift` is a per-axis half-grid offset
// (0 or 1); the double-mesh convention is
// q = (address * 2 + is_shift) / (mesh * 2).
namespace cppcrystal::kpoint {

// Linear grid-point index of a single-mesh address:
//   index = address[2]*mesh[0]*mesh[1] + address[1]*mesh[0] + address[0].
// The address must already be folded into [0, mesh).
[[nodiscard]] std::size_t grid_point_single_mesh(Vector3i const &address,
                                                 Vector3i const &mesh) noexcept;

// Linear grid-point index of a double-mesh address (halves the doubled address,
// folds modulo the mesh).
[[nodiscard]] std::size_t grid_point_double_mesh(Vector3i const &address_double,
                                                 Vector3i const &mesh) noexcept;

// Double-mesh address of a single-mesh address with a half-grid shift:
//   address_double = address * 2 + (is_shift != 0), folded to the
//   parallelepiped.
[[nodiscard]] Vector3i address_double_mesh(Vector3i const &address,
                                           Vector3i const &mesh,
                                           Vector3i const &is_shift) noexcept;

// All grid addresses of `mesh`, indexed by grid-point index (size
// mesh[0]*mesh[1]*mesh[2]).
[[nodiscard]] std::vector<Vector3i> all_grid_addresses(Vector3i const &mesh);

// Grid-point index of a grid address with no shift.
[[nodiscard]] std::size_t grid_point_from_address(Vector3i const &grid_address,
                                                  Vector3i const &mesh) noexcept;

} // namespace cppcrystal::kpoint
