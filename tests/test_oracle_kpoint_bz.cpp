// Oracle test for Brillouin-zone relocation (kpoint.c): relocate_BZ_grid_address
// and BZ_grid_points_by_rotations must reproduce
// spg_relocate_dense_BZ_grid_address / spg_get_dense_BZ_grid_points_by_rotations.

#include "oracle.hpp"

#include <spglib/kpoint/brillouin_zone.hpp>
#include <spglib/kpoint/grid.hpp>
#include <spglib/kpoint/reciprocal_mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace {

using spglib::Cell;
using spglib::Matrix3d;
using spglib::Matrix3i;
using spglib::Positions;
using spglib::Vector3d;
using spglib::Vector3i;

Matrix3d columns(Vector3d const &a, Vector3d const &b, Vector3d const &c) {
  Matrix3d l;
  l.col(0) = a;
  l.col(1) = b;
  l.col(2) = c;
  return l;
}

void check_relocate(Matrix3d const &real_lattice, Vector3i const &mesh,
                    Vector3i const &is_shift) {
  Matrix3d const rec = real_lattice.inverse().transpose();
  auto const grid = spglib::kpoint::all_grid_addresses(mesh);

  auto const got =
      spglib::kpoint::relocate_BZ_grid_address(grid, mesh, rec, is_shift);
  auto const ref = spglib::oracle::reference_relocate_BZ(grid, mesh, rec, is_shift);

  CHECK(got.bz_grid_address.size() == ref.num_bzgp);
  REQUIRE(got.bz_grid_address.size() == ref.bz_grid_address.size());
  for (std::size_t i = 0; i < ref.bz_grid_address.size(); ++i) {
    INFO("bz grid point " << i);
    CHECK(got.bz_grid_address[i] == ref.bz_grid_address[i]);
  }

  REQUIRE(got.bz_map.size() == ref.bz_map.size());
  for (std::size_t i = 0; i < ref.bz_map.size(); ++i) {
    INFO("bz_map index " << i);
    if (ref.bz_map[i] == ref.num_bzmesh) {
      CHECK_FALSE(got.bz_map[i].has_value());
    } else {
      REQUIRE(got.bz_map[i].has_value());
      CHECK(*got.bz_map[i] == ref.bz_map[i]);
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
    Matrix3d const fcc = columns(Vector3d(0, 2, 2), Vector3d(2, 0, 2),
                                 Vector3d(2, 2, 0));
    check_relocate(fcc, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
  SECTION("bcc primitive") {
    Matrix3d const bcc = columns(Vector3d(-2, 2, 2), Vector3d(2, -2, 2),
                                 Vector3d(2, 2, -2));
    check_relocate(bcc, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
  SECTION("hexagonal") {
    double const a = 3.0;
    double const c = 5.0;
    Matrix3d const hex = columns(Vector3d(a, 0, 0),
                                 Vector3d(-a / 2, a * std::sqrt(3.0) / 2, 0),
                                 Vector3d(0, 0, c));
    check_relocate(hex, Vector3i(4, 4, 4), Vector3i(0, 0, 0));
  }
}

TEST_CASE("BZ grid points by rotations match reference", "[oracle][kpoint]") {
  double const symprec = 1e-5;
  // A cubic cell supplies the reciprocal rotations and lattice.
  Positions p(1, 3);
  p.row(0) = Eigen::RowVector3d(0, 0, 0);
  Cell const cell(4.0 * Matrix3d::Identity(), p, {0});

  auto const ops = spglib::oracle::reference_symmetry(cell, symprec);
  std::vector<Matrix3i> rotations;
  for (auto const &op : ops)
    rotations.push_back(op.rotation);
  auto const rot_reciprocal =
      spglib::kpoint::point_group_reciprocal(rotations, false);

  Matrix3d const rec = cell.lattice().inverse().transpose();
  Vector3i const mesh(4, 4, 4);
  Vector3i const shift(0, 0, 0);
  auto const grid = spglib::kpoint::all_grid_addresses(mesh);
  auto const got_bz =
      spglib::kpoint::relocate_BZ_grid_address(grid, mesh, rec, shift);
  auto const ref_bz =
      spglib::oracle::reference_relocate_BZ(grid, mesh, rec, shift);

  for (Vector3i const &address : {Vector3i(1, 2, 0), Vector3i(2, 2, 0)}) {
    INFO("address " << address.transpose());
    auto const got = spglib::kpoint::BZ_grid_points_by_rotations(
        address, rot_reciprocal, mesh, shift, got_bz.bz_map);
    auto const ref = spglib::oracle::reference_BZ_grid_points_by_rotations(
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
