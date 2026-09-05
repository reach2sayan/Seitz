// Oracle test for the magnetic space-group determination (magnetic_spacegroup.c
// msg_identify_magnetic_space_group_type): the UNI number and MSG type found by
// magnetic::identify_magnetic_spacegroup_type must match the reference
// spg_get_magnetic_spacegroup_type_from_symmetry. The same magnetic symmetry
// (computed by the ported spin module) is fed to both sides.

#include "oracle.hpp"

#include "magnetic/identify.hpp"
#include "spin/search.hpp"
#include "symmetry/search.hpp"
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/operation_set.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <vector>

namespace {

using seitz::Cell;
using seitz::CollinearTensors;
using seitz::Lattice;
using seitz::MagneticCell;
using seitz::MagneticOperations;
using seitz::Matrix3d;
using seitz::noncollinear_tensors;
using seitz::Positions;
using seitz::SiteTensors;
using seitz::Vector3d;

Cell make_cell(double a, std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(Lattice{a * Matrix3d::Identity()}, p, types);
}

// The magnetic symmetry of `mcell`, via the ported spin module.
MagneticOperations magnetic_symmetry(MagneticCell const &input,
                                     bool with_time_reversal, bool is_axial,
                                     double symprec) {
  MagneticCell const mcell(input.cell(), input.tensors(),
                           is_axial ? seitz::TensorKind::axial
                                    : seitz::TensorKind::polar);
  seitz::symmetry::SymmetrySearch<seitz::GroupFamily::space> const
      spatial(mcell.cell(), {symprec});
  auto const sym_nonspin = spatial.operations();
  REQUIRE(sym_nonspin);
  seitz::spin::SpinSearch const spin_search(mcell, sym_nonspin.value(),
                                                 {{symprec}});
  auto const search =
      with_time_reversal
          ? spin_search.operations<seitz::TimeReversal::on>()
          : spin_search.operations<seitz::TimeReversal::off>();
  REQUIRE(search);
  return search->operations;
}

// Reference UNI number + MSG type via the magnetic symmetry directly.
std::pair<int, int> reference_uni(MagneticOperations const &ops,
                                  Matrix3d const &lattice, double symprec) {
  std::vector<int> rot(9 * ops.size());
  std::vector<double> trans(3 * ops.size());
  std::vector<int> timerev(ops.size());
  for (std::size_t s = 0; s < ops.size(); ++s) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        rot[9 * s + 3 * static_cast<std::size_t>(i) +
            static_cast<std::size_t>(j)] = ops[s].spatial.rotation(i, j);
      }
      trans[3 * s + static_cast<std::size_t>(i)] =
          ops[s].spatial.translation[i];
    }
    timerev[s] = ops[s].time_reversal ? 1 : 0;
  }
  double lat[3][3];
  seitz::oracle::to_c_lattice(lat, lattice);
  SpglibMagneticSpacegroupType const ref =
      spg_get_magnetic_spacegroup_type_from_symmetry(
          reinterpret_cast<int(*)[3][3]>(rot.data()),
          reinterpret_cast<double(*)[3]>(trans.data()), timerev.data(),
          static_cast<int>(ops.size()), lat, symprec);
  return {ref.uni_number, ref.type};
}

void check(MagneticCell const &mcell, bool with_time_reversal, bool is_axial,
           double symprec) {
  auto const ops =
      magnetic_symmetry(mcell, with_time_reversal, is_axial, symprec);
  seitz::magnetic::MagneticIdentification const identification(
      mcell.cell().lattice(), ops, {symprec});
  auto const got = identification.identify();
  REQUIRE(got);

  auto const [ref_uni, ref_type] =
      reference_uni(ops, mcell.cell().lattice().matrix(), symprec);
  REQUIRE(ref_uni != 0);
  REQUIRE(got->uni.value() == ref_uni);
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
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, -1.0}});
  SECTION("with time reversal") { check(mcell, true, true, 1e-5); }
  SECTION("without time reversal") { check(mcell, false, true, 1e-5); }
}

TEST_CASE("magnetic space-group type: non-collinear antiferromagnet",
          "[oracle][magnetic]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(
      cell, SiteTensors{noncollinear_tensors(
                {Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, -1.0)})});
  check(mcell, true, true, 1e-5);
}

TEST_CASE("magnetic space-group type: non-magnetic (type II)",
          "[oracle][magnetic]") {
  // All-equal spins with time reversal -> a grey (type-II) group.
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  MagneticCell const mcell(cell, SiteTensors{CollinearTensors{1.0, 1.0}});
  check(mcell, true, true, 1e-5);
}
