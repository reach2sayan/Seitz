// Oracle test for the k-point grid encoding (kgrid.c): grid_point_from_address
// must reproduce spg_get_dense_grid_point_from_address over every address in a
// window that exercises the boundary folding, for several meshes.

#include "oracle.hpp"

#include <cppcrystal/kpoint/grid.hpp>

#include <catch2/catch_test_macros.hpp>

namespace {

using cppcrystal::Vector3i;

void check_mesh(Vector3i const &mesh) {
  for (int x = -mesh[0]; x <= 2 * mesh[0]; ++x) {
    for (int y = -mesh[1]; y <= 2 * mesh[1]; ++y) {
      for (int z = -mesh[2]; z <= 2 * mesh[2]; ++z) {
        Vector3i const address(x, y, z);
        INFO("mesh " << mesh.transpose() << " address " << address.transpose());
        REQUIRE(cppcrystal::kpoint::grid_point_from_address(address, mesh) ==
                cppcrystal::oracle::reference_grid_point_from_address(address, mesh));
      }
    }
  }
}

} // namespace

TEST_CASE("grid_point_from_address matches reference", "[oracle][kpoint]") {
  check_mesh(Vector3i(4, 4, 4)); // even, isotropic
  check_mesh(Vector3i(3, 3, 3)); // odd
  check_mesh(Vector3i(2, 3, 4)); // anisotropic, mixed parity
}

TEST_CASE("all_grid_addresses indices are self-consistent", "[oracle][kpoint]") {
  for (Vector3i const &mesh :
       {Vector3i(4, 4, 4), Vector3i(3, 3, 3), Vector3i(2, 3, 4)}) {
    auto const grid = cppcrystal::kpoint::all_grid_addresses(mesh);
    REQUIRE(grid.size() ==
            static_cast<std::size_t>(mesh[0] * mesh[1] * mesh[2]));
    for (std::size_t gp = 0; gp < grid.size(); ++gp) {
      INFO("mesh " << mesh.transpose() << " grid point " << gp);
      REQUIRE(cppcrystal::kpoint::grid_point_from_address(grid[gp], mesh) == gp);
    }
  }
}
