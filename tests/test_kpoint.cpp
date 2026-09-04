// Unit tests for the k-point grid geometry and reciprocal point group (no
// reference spglib needed). Most of Mesh's own arithmetic is static_asserted in
// the header itself; what remains here is what needs a run.

#include <cppcrystal/kpoint/mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <vector>

using cppcrystal::Matrix3i;
using cppcrystal::kpoint::Address;
using cppcrystal::kpoint::Mesh;
using cppcrystal::kpoint::ReciprocalMesh;

TEST_CASE("Mesh index_of uses the a0-fastest encoding", "[kpoint]") {
  auto const mesh = Mesh::of({2, 3, 4}); // d0=2, d1=3, d0*d1=6
  REQUIRE(mesh);
  CHECK(mesh->index_of({0, 0, 0}) == 0U);
  CHECK(mesh->index_of({1, 0, 0}) == 1U); // a0 fastest
  CHECK(mesh->index_of({0, 1, 0}) == 2U); // +d0
  CHECK(mesh->index_of({0, 0, 1}) == 6U); // +d0*d1
  CHECK(mesh->index_of({1, 2, 3}) ==
        static_cast<std::size_t>(3 * 6 + 2 * 2 + 1));
}

TEST_CASE("Mesh addresses round-trip through index_of", "[kpoint]") {
  auto const mesh = Mesh::of({4, 4, 4});
  REQUIRE(mesh);
  REQUIRE(mesh->size() == 64U);
  std::size_t index = 0;
  for (Address const address : mesh->addresses()) {
    CHECK(mesh->index_of(address) == index);
    ++index;
  }
  CHECK(index == 64U);
}

TEST_CASE("a non-positive mesh has no Mesh", "[kpoint]") {
  CHECK_FALSE(Mesh::of({4, 0, 4}).has_value());
  CHECK_FALSE(Mesh::of({-1, 4, 4}).has_value());
  CHECK(Mesh::of({2, 2, 2})->size() == 8U);
}

TEST_CASE("reciprocal point group: transpose, inversion, dedup", "[kpoint]") {
  auto const mesh = Mesh::of({2, 2, 2});
  REQUIRE(mesh);
  std::vector<Matrix3i> const rots{Matrix3i::Identity()};

  auto const group = ReciprocalMesh::from_rotations(
      *mesh, rots, cppcrystal::TimeReversal::off);
  REQUIRE(group.rotations().size() == 1U);
  CHECK(group.rotations()[0] == Matrix3i::Identity());

  auto const grey =
      ReciprocalMesh::from_rotations(*mesh, rots, cppcrystal::TimeReversal::on);
  REQUIRE(grey.rotations().size() == 2U); // identity and its partner -I
  Matrix3i const minus_identity = -Matrix3i::Identity();
  CHECK((grey.rotations()[0] == minus_identity ||
         grey.rotations()[1] == minus_identity));
}
