// Oracle test for the irreducible / stabilized reciprocal mesh (kpoint.c):
// ir_reciprocal_mesh, stabilized_reciprocal_mesh, and grid_points_by_rotations
// must reproduce spg_get_dense_ir_reciprocal_mesh /
// spg_get_dense_stabilized_reciprocal_mesh / spg_get_dense_grid_points_by_
// rotations.

#include "oracle.hpp"

#include <cppcrystal/kpoint/reciprocal_mesh.hpp>
#include <cppcrystal/symmetry/pointgroup.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
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

Cell make_cell(Matrix3d const &lattice,
               std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i)
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  return Cell(Lattice{lattice}, p, types);
}

Matrix3d cubic(double a) { return a * Matrix3d::Identity(); }

Matrix3d hexagonal(double a, double c) {
  Matrix3d l;
  l.col(0) = Vector3d(a, 0.0, 0.0);
  l.col(1) = Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, c);
  return l;
}

void check_ir(Cell const &cell, Vector3i const &mesh, Vector3i const &is_shift,
              bool time_reversal, double symprec) {
  auto const got = cppcrystal::kpoint::ir_reciprocal_mesh(
      cell, mesh, is_shift,
      time_reversal ? cppcrystal::TimeReversal::on
                    : cppcrystal::TimeReversal::off,
      {symprec});
  REQUIRE(got);
  auto const ref = cppcrystal::oracle::reference_ir_reciprocal_mesh(
      cell, mesh, is_shift, time_reversal, symprec);

  CHECK(got->num_ir == ref.num_ir);
  REQUIRE(got->grid_address.size() == ref.grid_address.size());
  REQUIRE(got->ir_mapping_table.size() == ref.ir_mapping_table.size());
  for (std::size_t i = 0; i < ref.grid_address.size(); ++i) {
    INFO("grid point " << i);
    CHECK(got->grid_address[i] == ref.grid_address[i]);
    CHECK(got->ir_mapping_table[i] == ref.ir_mapping_table[i]);
  }
}

} // namespace

TEST_CASE("ir reciprocal mesh matches reference (normal path)",
          "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const pm3m = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  Cell const im3m =
      make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0});

  SECTION("gamma-centred") { check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), true, s); }
  SECTION("shifted") { check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(1, 1, 1), true, s); }
  SECTION("odd mesh") { check_ir(pm3m, Vector3i(3, 3, 3), Vector3i(0, 0, 0), true, s); }
  SECTION("no time reversal") { check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), false, s); }
  SECTION("bcc") { check_ir(im3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), true, s); }
}

TEST_CASE("ir reciprocal mesh matches reference (distortion path)",
          "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const p6mmm = make_cell(hexagonal(3.0, 5.0), {{0, 0, 0}}, {0});
  SECTION("gamma-centred") { check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 0), true, s); }
  SECTION("shifted c") { check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 1), true, s); }
  SECTION("no time reversal") { check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 0), false, s); }
}

TEST_CASE("stabilized reciprocal mesh matches reference", "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  auto const rotations = cppcrystal::symmetry::rotations_of(
      cppcrystal::oracle::reference_symmetry(cell, s));
  Vector3i const mesh(4, 4, 4);
  Vector3i const shift(0, 0, 0);

  auto compare = [&](std::vector<Vector3d> const &qpoints, bool tr) {
    auto const got = cppcrystal::kpoint::stabilized_reciprocal_mesh(
        mesh, shift,
        tr ? cppcrystal::TimeReversal::on : cppcrystal::TimeReversal::off,
        rotations, qpoints);
    auto const ref = cppcrystal::oracle::reference_stabilized_reciprocal_mesh(
        mesh, shift, tr, rotations, qpoints);
    CHECK(got.num_ir == ref.num_ir);
    REQUIRE(got.ir_mapping_table.size() == ref.ir_mapping_table.size());
    for (std::size_t i = 0; i < ref.ir_mapping_table.size(); ++i)
      CHECK(got.ir_mapping_table[i] == ref.ir_mapping_table[i]);
  };

  SECTION("gamma only (full symmetry)") {
    compare({Vector3d(0.0, 0.0, 0.0)}, true);
  }
  SECTION("q breaking symmetry") {
    compare({Vector3d(0.5, 0.0, 0.0)}, true);
  }
  SECTION("q set, no time reversal") {
    compare({Vector3d(0.5, 0.0, 0.0), Vector3d(0.0, 0.5, 0.0)}, false);
  }
}

TEST_CASE("grid points by rotations match reference", "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  auto const rotations = cppcrystal::symmetry::rotations_of(
      cppcrystal::oracle::reference_symmetry(cell, s));
  auto const rot_reciprocal =
      cppcrystal::kpoint::point_group_reciprocal(rotations,
                                             cppcrystal::TimeReversal::off);
  Vector3i const mesh(4, 4, 4);

  for (Vector3i const &shift : {Vector3i(0, 0, 0), Vector3i(1, 1, 1)}) {
    for (Vector3i const &address : {Vector3i(1, 2, 0), Vector3i(2, 0, 1)}) {
      INFO("shift " << shift.transpose() << " address " << address.transpose());
      auto const got = cppcrystal::kpoint::grid_points_by_rotations(
          address, rot_reciprocal, mesh, shift);
      auto const ref = cppcrystal::oracle::reference_grid_points_by_rotations(
          address, rot_reciprocal, mesh, shift);
      CHECK(got == ref);
    }
  }
}
