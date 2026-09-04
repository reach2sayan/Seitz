// Oracle test for Brillouin-zone relocation (kpoint.c):
// relocate_BZ_grid_address and BZ_grid_points_by_rotations must reproduce
// spg_relocate_dense_BZ_grid_address /
// spg_get_dense_BZ_grid_points_by_rotations.

#include "oracle.hpp"

#include <cppcrystal/kpoint/mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::Lattice;
using cppcrystal::Matrix3d;
using cppcrystal::Matrix3i;
using cppcrystal::Positions;
using cppcrystal::Vector3d;
using cppcrystal::Vector3i;

Matrix3d columns(Vector3d const &a, Vector3d const &b, Vector3d const &c) {
  Matrix3d l;
  l.col(0) = a;
  l.col(1) = b;
  l.col(2) = c;
  return l;
}

// The identity rotation is enough: BZ relocation itself does not use the
// reciprocal group, only its geometry.
cppcrystal::kpoint::ReciprocalMesh
reciprocal_of(Vector3i const &divisions, Vector3i const &is_shift,
              std::vector<Matrix3i> const &rots) {
  auto const mesh = cppcrystal::kpoint::Mesh::of(
      {divisions[0], divisions[1], divisions[2]},
      {is_shift[0] != 0, is_shift[1] != 0, is_shift[2] != 0});
  REQUIRE(mesh);
  return cppcrystal::kpoint::ReciprocalMesh::from_rotations(
      *mesh, rots, cppcrystal::TimeReversal::off);
}

std::vector<Vector3i> reference_grid(Vector3i const &divisions) {
  auto const mesh =
      cppcrystal::kpoint::Mesh::of({divisions[0], divisions[1], divisions[2]});
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
                       .brillouin_zone(cppcrystal::Lattice{rec});
  auto const ref = cppcrystal::oracle::reference_relocate_BZ(
      reference_grid(divisions), divisions, rec, is_shift);

  CHECK(got.addresses().size() == ref.num_bzgp);
  REQUIRE(got.addresses().size() == ref.bz_grid_address.size());
  for (std::size_t i = 0; i < ref.bz_grid_address.size(); ++i) {
    INFO("bz grid point " << i);
    auto const a = got.addresses()[i];
    CHECK(Vector3i(a[0], a[1], a[2]) == ref.bz_grid_address[i]);
  }

  for (std::size_t i = 0; i < ref.bz_map.size(); ++i) {
    INFO("bz_map index " << i);
    if (ref.bz_map[i] == ref.num_bzmesh) {
      CHECK_FALSE(got.map(i).has_value());
    } else {
      REQUIRE(got.map(i).has_value());
      CHECK(*got.map(i) == ref.bz_map[i]);
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
    double const a = 3.0;
    double const c = 5.0;
    Matrix3d const hex =
        columns(Vector3d(a, 0, 0), Vector3d(-a / 2, a * std::sqrt(3.0) / 2, 0),
                Vector3d(0, 0, c));
    check_relocate(hex, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
}

TEST_CASE("BZ grid points by rotations match reference", "[oracle][kpoint]") {
  double const symprec = 1e-5;
  // A cubic cell supplies the reciprocal rotations and lattice.
  Positions p(1, 3);
  p.row(0) = Eigen::RowVector3d(0, 0, 0);
  Cell const cell(Lattice{4.0 * Matrix3d::Identity()}, p, {0});

  auto const rotations =
      cppcrystal::oracle::reference_symmetry(cell, symprec).rotations();

  Matrix3d const rec = cell.lattice().matrix().inverse().transpose();
  Vector3i const mesh(4, 4, 4);
  Vector3i const shift(0, 0, 0);
  auto const reciprocal = reciprocal_of(mesh, shift, rotations);
  std::vector<Matrix3i> const rot_reciprocal(reciprocal.rotations().begin(),
                                             reciprocal.rotations().end());
  auto const got_bz = reciprocal.brillouin_zone(cppcrystal::Lattice{rec});
  auto const ref_bz = cppcrystal::oracle::reference_relocate_BZ(
      reference_grid(mesh), mesh, rec, shift);

  for (Vector3i const &address : {Vector3i(1, 2, 0), Vector3i(2, 2, 0)}) {
    INFO("address " << address.transpose());
    auto const got = got_bz.images_of({address[0], address[1], address[2]});
    auto const ref = cppcrystal::oracle::reference_BZ_grid_points_by_rotations(
        address, rot_reciprocal, mesh, shift, ref_bz.bz_map);

    REQUIRE(got.size() == ref.size());
    for (std::size_t i = 0; i < ref.size(); ++i) {
      if (ref[i] == ref_bz.num_bzmesh) {
        CHECK_FALSE(got[i].has_value());
      } else {
        REQUIRE(got[i].has_value());
        CHECK(*got[i] == ref[i]);
      }
    }
  }
}
