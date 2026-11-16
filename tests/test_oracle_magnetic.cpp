// Oracle test for the magnetic space-group determination (magnetic_spacegroup.c
// msg_identify_magnetic_space_group_type): the UNI number and MSG type found by
// MagneticOperations::spacegroup must match the reference
// spg_get_magnetic_spacegroup_type_from_symmetry. The same magnetic symmetry
// (computed by the ported spin module) is fed to both sides.

#include "oracle.hpp"
#include "structures.hpp"

#include "spin/search.hpp"
#include "symmetry/search.hpp"
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/operation_set.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <random>
#include <ranges>
#include <vector>

namespace {

using seitz::test::make_cell;

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

// The magnetic symmetry of `mcell`, via the ported spin module.
MagneticOperations magnetic_symmetry(MagneticCell const &input,
                                     bool with_time_reversal, bool is_axial,
                                     double symprec) {
  MagneticCell const mcell(input.cell(), input.tensors(),
                           is_axial ? seitz::TensorKind::axial
                                    : seitz::TensorKind::polar);
  seitz::symmetry::SymmetrySearch<seitz::GroupFamily::space> const spatial(
      mcell.cell(), {symprec});
  auto const sym_nonspin = spatial.operations();
  REQUIRE(sym_nonspin);
  seitz::spin::SpinSearch const spin_search(mcell, sym_nonspin.value(),
                                            {{symprec}});
  auto const search = with_time_reversal
                          ? spin_search.operations<seitz::TimeReversal::on>()
                          : spin_search.operations<seitz::TimeReversal::off>();
  REQUIRE(search);
  return search->operations;
}

// Reference UNI number + MSG type via the magnetic symmetry directly.
std::pair<int, int> reference_uni(MagneticOperations const &ops,
                                  Matrix3d const &lattice, double symprec) {
  seitz::oracle::COperations c_ops(
      ops | std::views::transform(&seitz::MagneticSymmetryOperation::spatial));
  std::vector<int> timerev{
      std::from_range,
      ops | std::views::transform([](auto const &op) -> int {
        return op.time_reversal ? 1 : 0;
      })};
  double lat[3][3];
  seitz::oracle::to_c_lattice(lat, lattice);
  SpglibMagneticSpacegroupType const ref =
      spg_get_magnetic_spacegroup_type_from_symmetry(
          c_ops.rotations(), c_ops.translations(), timerev.data(),
          c_ops.size(), lat, symprec);
  return {ref.uni_number, ref.type};
}

void check(MagneticCell const &mcell, bool with_time_reversal, bool is_axial,
           double symprec) {
  auto const ops =
      magnetic_symmetry(mcell, with_time_reversal, is_axial, symprec);
  auto const got = ops.spacegroup(mcell.cell().lattice(), {symprec});
  REQUIRE(got);

  auto const [ref_uni, ref_type] =
      reference_uni(ops, mcell.cell().lattice().matrix(), symprec);
  REQUIRE(ref_uni != 0);
  REQUIRE(got->uni.value() == ref_uni);
  REQUIRE(static_cast<int>(got->type) == ref_type);
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

TEST_CASE("magnetic space-group type: a jittered antiferromagnetic supercell",
          "[oracle][magnetic]") {
  // A 2x2x2 supercell whose coordinates are displaced well below symprec, so
  // the search still resolves the full symmetry but the operations carry
  // noise. Transforming the eight pure translations into the reference setting
  // multiplies that noise by det(tmat) = 8 and lands some of them just below
  // zero, i.e. just below 1 after wrapping: deduplicating them componentwise
  // rather than modulo 1 counted one translation as four and failed the
  // determination outright.
  constexpr double kSymprec = 1e-3;
  constexpr double kJitter = 1e-4; // fractional, two decades below symprec
  constexpr int k = 2;

  std::vector<std::array<double, 3>> positions;
  std::vector<int> types;
  CollinearTensors moments;
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> jitter(-kJitter, kJitter);
  for (auto const [a, b, c, site] : std::views::cartesian_product(
           std::views::iota(0, k), std::views::iota(0, k),
           std::views::iota(0, k), std::views::iota(0, 2))) {
    double const offset = 0.5 * site;
    positions.push_back({(a + offset) / k + jitter(rng),
                         (b + offset) / k + jitter(rng),
                         (c + offset) / k + jitter(rng)});
    types.push_back(site);
    moments.push_back(site == 0 ? 1.0 : -1.0);
  }

  Cell const cell = make_cell(4.0 * k, positions, types);
  MagneticCell const mcell(cell, SiteTensors{moments});
  check(mcell, true, true, kSymprec);
}
