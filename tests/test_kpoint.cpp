// Unit tests for the k-point grid encoding and reciprocal point group (no
// reference spglib needed).

#include <spglib/kpoint/grid.hpp>
#include <spglib/kpoint/reciprocal_mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <vector>

using spglib::Matrix3i;
using spglib::Vector3i;

TEST_CASE("grid_point_single_mesh uses the default encoding", "[kpoint]") {
  Vector3i const mesh(2, 3, 4); // m0=2, m1=3, m0*m1=6
  using spglib::kpoint::grid_point_single_mesh;
  REQUIRE(grid_point_single_mesh(Vector3i(0, 0, 0), mesh) == 0U);
  REQUIRE(grid_point_single_mesh(Vector3i(1, 0, 0), mesh) == 1U); // a0 fastest
  REQUIRE(grid_point_single_mesh(Vector3i(0, 1, 0), mesh) == 2U); // +m0
  REQUIRE(grid_point_single_mesh(Vector3i(0, 0, 1), mesh) == 6U); // +m0*m1
  REQUIRE(grid_point_single_mesh(Vector3i(1, 2, 3), mesh) ==
          static_cast<std::size_t>(3 * 6 + 2 * 2 + 1));
}

TEST_CASE("all_grid_addresses round-trips through grid_point_from_address",
          "[kpoint]") {
  Vector3i const mesh(4, 4, 4);
  auto const grid = spglib::kpoint::all_grid_addresses(mesh);
  REQUIRE(grid.size() == 64U);
  for (std::size_t gp = 0; gp < grid.size(); ++gp) {
    REQUIRE(spglib::kpoint::grid_point_from_address(grid[gp], mesh) == gp);
  }
}

TEST_CASE("reciprocal point group: transpose, inversion, dedup", "[kpoint]") {
  std::vector<Matrix3i> const rots{Matrix3i::Identity()};
  auto const group = spglib::kpoint::point_group_reciprocal(rots, false);
  REQUIRE(group.size() == 1U);
  REQUIRE(group[0] == Matrix3i::Identity());

  auto const grey = spglib::kpoint::point_group_reciprocal(rots, true);
  REQUIRE(grey.size() == 2U); // identity and its inversion partner -I
  Matrix3i const minus_identity = -Matrix3i::Identity();
  REQUIRE((grey[0] == minus_identity || grey[1] == minus_identity));
}
