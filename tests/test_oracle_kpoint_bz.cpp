// Oracle test for Brillouin-zone relocation (kpoint.c):
// relocate_BZ_grid_address and BZ_grid_points_by_rotations must reproduce
// spg_relocate_dense_BZ_grid_address /
// spg_get_dense_BZ_grid_points_by_rotations.

#include "oracle.hpp"
#include "structures.hpp"

#include <seitz/kpoint/mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace {

using seitz::test::columns;

using seitz::Cell;
using seitz::Lattice;
using seitz::Matrix3d;
using seitz::Matrix3i;
using seitz::Positions;
using seitz::Vector3d;
using seitz::Vector3i;

// The identity rotation is enough: BZ relocation itself does not use the
// reciprocal group, only its geometry.
seitz::kpoint::ReciprocalMesh
reciprocal_of(Vector3i const &divisions, Vector3i const &is_shift,
              std::vector<Matrix3i> const &rots) {
  auto const mesh = seitz::kpoint::Mesh::of(
      {divisions[0], divisions[1], divisions[2]},
      {is_shift[0] != 0, is_shift[1] != 0, is_shift[2] != 0});
  REQUIRE(mesh);
  return seitz::kpoint::ReciprocalMesh::from_rotations(
      *mesh, rots, seitz::TimeReversal::off);
}

std::vector<Vector3i> reference_grid(Vector3i const &divisions) {
  auto const mesh =
      seitz::kpoint::Mesh::of({divisions[0], divisions[1], divisions[2]});
  std::vector<Vector3i> out;
  for (auto const a : mesh->addresses())
    out.emplace_back(a[0], a[1], a[2]);
  return out;
}

void check_relocate(Matrix3d const &real_lattice, Vector3i const &divisions,
                    Vector3i const &is_shift) {
  Matrix3d const rec = real_lattice.inverse().transpose();
  std::vector<Matrix3i> const identity{Matrix3i::Identity()};
  auto const got = reciprocal_of(divisions, is_shift, identity)
                       .brillouin_zone(seitz::Lattice{rec});
  auto const ref = seitz::oracle::reference_relocate_BZ(
      reference_grid(divisions), divisions, rec, is_shift);

  CHECK(got.addresses().size() == ref.num_bzgp);
  REQUIRE(got.addresses().size() == ref.bz_grid_address.size());
  for (auto const [i, point] :
       std::views::zip(got.addresses(), ref.bz_grid_address) |
           std::views::enumerate) {
    auto const &[a, ref_address] = point;
    INFO("bz grid point " << static_cast<std::size_t>(i));
    CHECK(Vector3i(a[0], a[1], a[2]) == ref_address);
  }

  for (auto const [i, ref_mapped] : ref.bz_map | std::views::enumerate) {
    auto const index = static_cast<std::size_t>(i);
    INFO("bz_map index " << index);
    if (ref_mapped == ref.num_bzmesh) {
      CHECK_FALSE(got.map(index).has_value());
    } else {
      REQUIRE(got.map(index).has_value());
      CHECK(*got.map(index) == ref_mapped);
    }
  }
}

} // namespace

TEST_CASE("BZ relocation matches reference", "[oracle][kpoint]") {
  SECTION("simple cubic") {
    check_relocate(4.0 * Matrix3d::Identity(), Vector3i(4, 4, 4),
                   Vector3i(0, 0, 0));
  }
  SECTION("simple cubic, shifted") {
    check_relocate(4.0 * Matrix3d::Identity(), Vector3i(4, 4, 4),
                   Vector3i(1, 1, 1));
  }
  SECTION("fcc primitive") {
    Matrix3d const fcc =
        columns(Vector3d(0, 2, 2), Vector3d(2, 0, 2), Vector3d(2, 2, 0));
    check_relocate(fcc, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
  SECTION("bcc primitive") {
    Matrix3d const bcc =
        columns(Vector3d(-2, 2, 2), Vector3d(2, -2, 2), Vector3d(2, 2, -2));
    check_relocate(bcc, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
  SECTION("hexagonal") {
    check_relocate(seitz::test::hexagonal(3.0, 5.0), Vector3i(4, 4, 4),
                   Vector3i(0, 0, 0));
  }
}

TEST_CASE("BZ grid points by rotations match reference", "[oracle][kpoint]") {
  double const symprec = 1e-5;
  // A cubic cell supplies the reciprocal rotations and lattice.
  Positions p(1, 3);
  p.row(0) = Eigen::RowVector3d(0, 0, 0);
  Cell const cell(Lattice{4.0 * Matrix3d::Identity()}, p, {0});

  auto const rotations =
      seitz::oracle::reference_symmetry(cell, symprec).rotations();

  Matrix3d const rec = cell.lattice().matrix().inverse().transpose();
  Vector3i const mesh(4, 4, 4);
  Vector3i const shift(0, 0, 0);
  auto const reciprocal = reciprocal_of(mesh, shift, rotations);
  std::vector<Matrix3i> const rot_reciprocal(reciprocal.rotations().begin(),
                                             reciprocal.rotations().end());
  auto const got_bz = reciprocal.brillouin_zone(seitz::Lattice{rec});
  auto const ref_bz = seitz::oracle::reference_relocate_BZ(
      reference_grid(mesh), mesh, rec, shift);

  for (Vector3i const &address : {Vector3i(1, 2, 0), Vector3i(2, 2, 0)}) {
    INFO("address " << address.transpose());
    auto const got = got_bz.images_of({address[0], address[1], address[2]});
    auto const ref = seitz::oracle::reference_BZ_grid_points_by_rotations(
        address, rot_reciprocal, mesh, shift, ref_bz.bz_map);

    REQUIRE(got.size() == ref.size());
    for (auto const &[image, ref_image] : std::views::zip(got, ref)) {
      if (ref_image == ref_bz.num_bzmesh) {
        CHECK_FALSE(image.has_value());
      } else {
        REQUIRE(image.has_value());
        CHECK(*image == ref_image);
      }
    }
  }
}
