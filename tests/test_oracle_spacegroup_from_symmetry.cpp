// Oracle test for the symmetry -> space-group pipeline (the foundation the
// magnetic determination builds on): OperationSet::spacegroup
// (= prm_get_primitive_symmetry + spa_search_spacegroup_with_symmetry) must
// reproduce spg_get_spacegroup_type_from_symmetry. The same operation set is
// fed to both sides, so this tests the pipeline, not the symmetry search.

#include "oracle.hpp"
#include <cppcrystal/core/operation_set.hpp>

#include "spacegroup/spacegroup.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::Lattice;
using cppcrystal::Matrix3d;
using cppcrystal::Operations;
using cppcrystal::Positions;
using cppcrystal::oracle::reference_symmetry;

Cell make_cell(Matrix3d const &lattice,
               std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(Lattice{lattice}, p, types);
}

Matrix3d cubic(double a) { return a * Matrix3d::Identity(); }

Matrix3d hexagonal(double a, double c) {
  Matrix3d l;
  l.col(0) = Eigen::Vector3d(a, 0.0, 0.0);
  l.col(1) = Eigen::Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0);
  l.col(2) = Eigen::Vector3d(0.0, 0.0, c);
  return l;
}

Matrix3d tetragonal(double a, double c) {
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = a;
  l(1, 1) = a;
  l(2, 2) = c;
  return l;
}

void to_c_operations(std::vector<int> &rot, std::vector<double> &trans,
                     Operations const &ops) {
  rot.resize(9 * ops.size());
  trans.resize(3 * ops.size());
  for (std::size_t s = 0; s < ops.size(); ++s) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rot[9 * s + 3 * static_cast<std::size_t>(i) +
            static_cast<std::size_t>(j)] = ops[s].rotation(i, j);
      }
      trans[3 * s + static_cast<std::size_t>(i)] = ops[s].translation[i];
    }
  }
}

// Compare port OperationSet::spacegroup<conventional> against the reference
// spg_get_spacegroup_type_from_symmetry, both fed the same operation set.
void check_conventional(Cell const &cell, double symprec) {
  auto const ops = reference_symmetry(cell, symprec);
  REQUIRE(!ops.empty());

  auto const got = ops.spacegroup(cell.lattice().matrix(), {symprec});
  REQUIRE(got);

  std::vector<int> rot;
  std::vector<double> trans;
  to_c_operations(rot, trans, ops);
  double lat[3][3];
  cppcrystal::oracle::to_c_lattice(lat, cell.lattice().matrix());
  SpglibSpacegroupType const ref = spg_get_spacegroup_type_from_symmetry(
      reinterpret_cast<int(*)[3][3]>(rot.data()),
      reinterpret_cast<double(*)[3]>(trans.data()),
      static_cast<int>(ops.size()), lat, symprec);

  REQUIRE(ref.number != 0);
  REQUIRE(got->type().number == ref.number);
  REQUIRE(got->hall.index() == ref.hall_number);
  REQUIRE(got->type().international_short ==
          std::string(ref.international_short));
}

} // namespace

TEST_CASE("OperationSet::spacegroup matches reference (conventional)",
          "[oracle][spacegroup][from_symmetry]") {
  double const s = 1e-5;
  SECTION("Pm-3m (primitive cubic)") {
    check_conventional(make_cell(cubic(4.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("Im-3m (body-centred cubic)") {
    check_conventional(
        make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0}), s);
  }
  SECTION("Fm-3m (face-centred cubic)") {
    check_conventional(
        make_cell(cubic(4.0),
                  {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}},
                  {0, 0, 0, 0}),
        s);
  }
  SECTION("P6/mmm (hexagonal)") {
    check_conventional(make_cell(hexagonal(3.0, 5.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("P4/mmm (tetragonal)") {
    check_conventional(make_cell(tetragonal(3.0, 5.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("rutile-like (tetragonal, 2 species)") {
    check_conventional(make_cell(tetragonal(4.6, 3.0),
                                 {{0, 0, 0},
                                  {0.5, 0.5, 0.5},
                                  {0.3, 0.3, 0.0},
                                  {0.7, 0.7, 0.0},
                                  {0.8, 0.2, 0.5},
                                  {0.2, 0.8, 0.5}},
                                 {0, 0, 1, 1, 1, 1}),
                       s);
  }
}

// Exercises the LatticeSetting::primitive branch (identity primitive lattice),
// mirroring spg_get_hall_number_from_symmetry.
TEST_CASE("OperationSet::spacegroup matches reference (primitive)",
          "[oracle][spacegroup][from_symmetry]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0});
  auto const ops = reference_symmetry(cell, s);
  REQUIRE(!ops.empty());

  std::vector<int> rot;
  std::vector<double> trans;
  to_c_operations(rot, trans, ops);
  int const ref_hall = spg_get_hall_number_from_symmetry(
      reinterpret_cast<int(*)[3][3]>(rot.data()),
      reinterpret_cast<double(*)[3]>(trans.data()),
      static_cast<int>(ops.size()), s);

  auto const got = ops.spacegroup<cppcrystal::LatticeSetting::primitive>(
      Matrix3d::Identity(), {s});

  REQUIRE((ref_hall == 0) == (!got.has_value()));
  if (got) {
    REQUIRE(got->hall.index() == ref_hall);
  }
}
