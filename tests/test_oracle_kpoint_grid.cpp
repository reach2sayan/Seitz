// Oracle test for the k-point grid encoding (kgrid.c): Mesh::index_of must
// reproduce spg_get_dense_grid_point_from_address over every address in a
// window that exercises the boundary folding, for several meshes.

#include "oracle.hpp"

#include <seitz/kpoint/mesh.hpp>

#include <catch2/catch_test_macros.hpp>

namespace {

using seitz::Vector3i;
using seitz::kpoint::Address;
using seitz::kpoint::Mesh;

void check_mesh(Address const &divisions) {
  auto const mesh = Mesh::of(divisions);
  REQUIRE(mesh);
  Vector3i const as_vector(divisions[0], divisions[1], divisions[2]);
  for (int x = -divisions[0]; x <= 2 * divisions[0]; ++x) {
    for (int y = -divisions[1]; y <= 2 * divisions[1]; ++y) {
      for (int z = -divisions[2]; z <= 2 * divisions[2]; ++z) {
        INFO("mesh " << as_vector.transpose() << " address " << x << " " << y
                     << " " << z);
        REQUIRE(mesh->index_of({x, y, z}) ==
                seitz::oracle::reference_grid_point_from_address(
                    Vector3i(x, y, z), as_vector));
      }
    }
  }
}

} // namespace

TEST_CASE("Mesh::index_of matches the reference encoding", "[oracle][kpoint]") {
  check_mesh({4, 4, 4}); // even, isotropic
  check_mesh({3, 3, 3}); // odd
  check_mesh({2, 3, 4}); // anisotropic, mixed parity
}

TEST_CASE("Mesh addresses are self-consistent", "[oracle][kpoint]") {
  for (Address const &divisions :
       {Address{4, 4, 4}, Address{3, 3, 3}, Address{2, 3, 4}}) {
    auto const mesh = Mesh::of(divisions);
    REQUIRE(mesh);
    REQUIRE(mesh->size() == static_cast<std::size_t>(
                                divisions[0] * divisions[1] * divisions[2]));
    std::size_t index = 0;
    for (Address const address : mesh->addresses()) {
      INFO("grid point " << index);
      REQUIRE(mesh->index_of(address) == index);
      ++index;
    }
  }
}
