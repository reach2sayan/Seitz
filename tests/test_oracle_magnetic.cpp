// Oracle test for the magnetic space-group determination (magnetic_spacegroup.c
// msg_identify_magnetic_space_group_type): the UNI number and MSG type found by
// magnetic::identify_magnetic_spacegroup_type must match the reference
// spg_get_magnetic_spacegroup_type_from_symmetry. The same magnetic symmetry
// (computed by the ported spin module) is fed to both sides.

#include "oracle.hpp"

#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/magnetic/magnetic_spacegroup.hpp>
#include <cppcrystal/spin/spin.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::CollinearTensors;
using cppcrystal::MagneticCell;
using cppcrystal::MagneticSymmetryOperations;
using cppcrystal::Matrix3d;
using cppcrystal::noncollinear_tensors;
using cppcrystal::Positions;
using cppcrystal::SiteTensors;
using cppcrystal::Vector3d;

Cell make_cell(double a, std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(a * Matrix3d::Identity(), p, types);
}

// The magnetic symmetry of `mcell`, via the ported spin module.
MagneticSymmetryOperations magnetic_symmetry(MagneticCell const &mcell,
                                             bool with_time_reversal,
                                             bool is_axial, double symprec) {
  auto const sym_nonspin =
      cppcrystal::symmetry::find_symmetry(mcell.cell(), symprec);
  REQUIRE(sym_nonspin);
  auto const search = cppcrystal::spin::operations_with_site_tensors(
      sym_nonspin.value(), mcell, with_time_reversal, is_axial, symprec);
  REQUIRE(search);
  return search->operations;
}

// Reference UNI number + MSG type via the magnetic symmetry directly.
std::pair<int, int>
reference_uni(MagneticSymmetryOperations const &ops, Matrix3d const &lattice,
              double symprec) {
  std::vector<int> rot(9 * ops.size());
  std::vector<double> trans(3 * ops.size());
  std::vector<int> timerev(ops.size());
  for (std::size_t s = 0; s < ops.size(); ++s) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rot[9 * s + 3 * static_cast<std::size_t>(i) +
            static_cast<std::size_t>(j)] = ops[s].rotation(i, j);
      }
      trans[3 * s + static_cast<std::size_t>(i)] = ops[s].translation[i];
    }
    timerev[s] = ops[s].time_reversal ? 1 : 0;
  }
  double lat[3][3];
  cppcrystal::oracle::to_c_lattice(lat, lattice);
  SpglibMagneticSpacegroupType const ref =
      spg_get_magnetic_spacegroup_type_from_symmetry(
          reinterpret_cast<int(*)[3][3]>(rot.data()),
          reinterpret_cast<double(*)[3]>(trans.data()), timerev.data(),
          static_cast<int>(ops.size()), lat, symprec);
  return {ref.uni_number, ref.type};
}

void check(MagneticCell const &mcell, bool with_time_reversal, bool is_axial,
           double symprec) {
  auto const ops = magnetic_symmetry(mcell, with_time_reversal, is_axial, symprec);
  auto const got = cppcrystal::magnetic::identify_magnetic_spacegroup_type(
      mcell.cell().lattice(), ops, symprec);
  REQUIRE(got);

  auto const [ref_uni, ref_type] =
      reference_uni(ops, mcell.cell().lattice(), symprec);
  REQUIRE(ref_uni != 0);
  REQUIRE(got->uni_number == ref_uni);
  REQUIRE(static_cast<int>(got->msg_type) == ref_type);
}

} // namespace

TEST_CASE("magnetic space-group type: collinear ferromagnet",
          "[oracle][magnetic]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}}, {0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0}});
  check(mcell, true, true, 1e-5);
}

TEST_CASE("magnetic space-group type: collinear antiferromagnet",
          "[oracle][magnetic]") {
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, -1.0}});
  SECTION("with time reversal") { check(mcell, true, true, 1e-5); }
  SECTION("without time reversal") { check(mcell, false, true, 1e-5); }
}

TEST_CASE("magnetic space-group type: non-collinear antiferromagnet",
          "[oracle][magnetic]") {
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(
      cell, SiteTensors{noncollinear_tensors(
                {Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, -1.0)})});
  check(mcell, true, true, 1e-5);
}

TEST_CASE("magnetic space-group type: non-magnetic (type II)",
          "[oracle][magnetic]") {
  // All-equal spins with time reversal -> a grey (type-II) group.
  Cell const cell =
      make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, 1.0}});
  check(mcell, true, true, 1e-5);
}
