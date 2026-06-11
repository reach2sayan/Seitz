#include <spglib/kpoint/grid.hpp>

#include <cstddef>

namespace spglib::kpoint {

namespace {

// A mesh must be strictly positive on every axis; otherwise the modulo / index
// arithmetic below is undefined (mod-by-zero) or allocates wildly. Guarded at
// the two roots (all_grid_addresses, grid_point_double_mesh) so every grid /
// reciprocal-mesh function degrades to empty/zero instead of UB on a bad mesh.
[[nodiscard]] bool valid_mesh(Vector3i const &mesh) noexcept {
  return mesh[0] > 0 && mesh[1] > 0 && mesh[2] > 0;
}

// Euclidean modulo into [0, m): kgrid.c modulo_i3.
void modulo_i3(Vector3i &v, Vector3i const &m) noexcept {
  for (int i = 0; i < 3; ++i) {
    v[i] %= m[i];
    if (v[i] < 0) {
      v[i] += m[i];
    }
  }
}

// Fold a single-mesh address onto the parallelepiped (default branch, no
// GRID_BOUNDARY_AS_NEGATIVE): kgrid.c reduce_grid_address.
void inline reduce_grid_address(Vector3i &address, Vector3i const &mesh) noexcept {
  address.array() -=
      mesh.array() * (address.array() > mesh.array() / 2).cast<int>();
}

// Fold a double-mesh address onto the parallelepiped: reduce_grid_address_double.
void inline reduce_grid_address_double(Vector3i &address,
                                Vector3i const &mesh) noexcept {
  address.array() -=
      2 * mesh.array() * (address.array() > mesh.array()).cast<int>();
}

} // namespace

std::size_t grid_point_single_mesh(Vector3i const &address,
                                   Vector3i const &mesh) noexcept {
  // The address is folded into [0, mesh) before this is called, so the casts to
  // size_t are non-negative (and avoid int overflow on large meshes).
  const std::size_t sx = 1;
  const std::size_t sy = static_cast<std::size_t>(mesh[0]);
  const std::size_t sz =
      static_cast<std::size_t>(mesh[0]) * static_cast<std::size_t>(mesh[1]);

  return static_cast<std::size_t>(address[0]) * sx +
         static_cast<std::size_t>(address[1]) * sy +
         static_cast<std::size_t>(address[2]) * sz;
}

std::size_t grid_point_double_mesh(Vector3i const &address_double,
                                   Vector3i const &mesh) noexcept {
  if (!valid_mesh(mesh)) {
    return 0;
  }
  Vector3i address =
      (address_double.cast<double>().array() / 2.0).floor().cast<int>();
  modulo_i3(address, mesh);
  return grid_point_single_mesh(address, mesh);
}

Vector3i address_double_mesh(Vector3i const &address, Vector3i const &mesh,
                             Vector3i const &is_shift) noexcept {
  Vector3i address_double =
      2 * address.array() + (is_shift.array() != 0).cast<int>();
  reduce_grid_address_double(address_double, mesh);
  return address_double;
}

std::vector<Vector3i> all_grid_addresses(Vector3i const &mesh) {
  if (!valid_mesh(mesh)) {
    return {};
  }
  const auto n = static_cast<std::size_t>(mesh.prod());
  std::vector<Vector3i> grid_address(n);
  for (int i = 0; i < mesh[0]; ++i) {
    for (int j = 0; j < mesh[1]; ++j) {
      for (int k = 0; k < mesh[2]; ++k) {
        Vector3i const address(i, j, k);
        Vector3i reduced = address;
        reduce_grid_address(reduced, mesh);
        grid_address[grid_point_single_mesh(address, mesh)] = reduced;
      }
    }
  }
  return grid_address;
}

std::size_t grid_point_from_address(Vector3i const &grid_address,
                                    Vector3i const &mesh) noexcept {
  Vector3i const address_double =
      address_double_mesh(grid_address, mesh, Vector3i::Zero());
  return grid_point_double_mesh(address_double, mesh);
}

} // namespace spglib::kpoint
