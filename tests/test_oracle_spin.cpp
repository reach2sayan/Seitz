// Oracle test for the magnetic symmetry search (Phase 7, spin.c): the magnetic
// operations / equivalent atoms / primitive lattice found by
// spin::operations_with_site_tensors must match the reference
// spg_get_symmetry_with_site_tensors.

#include "oracle.hpp"
#include "structures.hpp"

#include "spin/search.hpp"
#include "symmetry/search.hpp"
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/operation_set.hpp>

#include "math/integer_matrix.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using seitz::test::make_cell;

using seitz::Cell;
using seitz::CollinearTensors;
using seitz::Lattice;
using seitz::MagneticCell;
using seitz::MagneticOperations;
using seitz::Matrix3d;
using seitz::Matrix3i;
using seitz::noncollinear_tensors;
using seitz::Positions;
using seitz::SiteTensors;
using seitz::Vector3d;

// Reference magnetic symmetry via spg_get_symmetry_with_site_tensors. `tensors`
// is flat: collinear (rank 0) = one per atom, non-collinear (rank 1) = three
// per atom. Returns operations with time_reversal derived from spin_flips
// (spin_flips = 1 - 2*timerev, so timerev <=> spin_flips == -1).
struct Reference {
  std::vector<seitz::MagneticSymmetryOperation> operations;
  std::vector<int> equivalent_atoms;
  Matrix3d primitive_lattice;
};

Reference reference_magnetic(MagneticCell const &mcell,
                             bool with_time_reversal, bool is_axial,
                             double symprec) {
  seitz::oracle::CCell c(mcell.cell());
  std::vector<double> const tensors = seitz::oracle::flat_tensors(mcell);
  int const tensor_rank = seitz::oracle::tensor_rank(mcell);
  int const n = c.num_atom();
  int const max_size = 384;
  seitz::oracle::COperations buffer(static_cast<std::size_t>(max_size));
  std::vector<int> equiv(static_cast<std::size_t>(n));
  std::vector<int> spin_flips(static_cast<std::size_t>(max_size));
  double prim[3][3];

  int const num = spg_get_symmetry_with_site_tensors(
      buffer.rotations(), buffer.translations(), equiv.data(), prim,
      spin_flips.data(), max_size, c.lattice, c.pos(), c.types.data(),
      tensors.data(), tensor_rank, n, with_time_reversal, is_axial, symprec);
  REQUIRE(num > 0);

  Reference ref;
  auto const spatial = buffer.operations(num);
  ref.operations = {
      std::from_range,
      std::views::zip(spatial, spin_flips) |
          std::views::transform([](auto const &entry) {
            auto const &[op, flip] = entry;
            return seitz::MagneticSymmetryOperation{op, flip == -1};
          })};
  ref.equivalent_atoms = equiv;
  seitz::oracle::from_c_lattice(ref.primitive_lattice, prim);
  return ref;
}


void check(MagneticCell const &input, bool with_time_reversal, bool is_axial,
           double symprec) {
  MagneticCell const mcell(input.cell(), input.tensors(),
                           is_axial ? seitz::TensorKind::axial
                                    : seitz::TensorKind::polar);
  seitz::symmetry::SymmetrySearch<seitz::GroupFamily::space> const
      search(mcell.cell(), {symprec});
  auto const sym_nonspin = search.operations();
  REQUIRE(sym_nonspin);

  seitz::spin::SpinSearch const spin_search(mcell, sym_nonspin.value(),
                                                 {{symprec}});
  auto const got =
      with_time_reversal
          ? spin_search.operations<seitz::TimeReversal::on>()
          : spin_search.operations<seitz::TimeReversal::off>();
  REQUIRE(got);

  auto const ref =
      reference_magnetic(mcell, with_time_reversal, is_axial, symprec);

  // Operations match as a set (the spatial-symmetry order may differ).
  REQUIRE(got->operations.size() == ref.operations.size());
  for (auto const &op : ref.operations) {
    INFO("missing reference operation");
    REQUIRE(std::ranges::any_of(got->operations, [&](auto const &g) {
      return g.time_reversal == op.time_reversal &&
             seitz::same_operation(g.spatial, op.spatial, symprec);
    }));
  }
  REQUIRE(got->equivalent_atoms == ref.equivalent_atoms);
  REQUIRE(seitz::test::same_lattice(got->primitive_lattice,
                                    ref.primitive_lattice));
}

} // namespace

TEST_CASE("collinear ferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}}, {0});
  CollinearTensors const spins{1.0};
  MagneticCell const mcell(cell, SiteTensors{spins});
  check(mcell, true, true, 1e-5);
}

TEST_CASE("collinear antiferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  CollinearTensors const spins{1.0, -1.0};
  MagneticCell const mcell(cell, SiteTensors{spins});
  SECTION("with time reversal (family group)") {
    check(mcell, true, true, 1e-5);
  }
  SECTION("without time reversal (maximal subgroup)") {
    check(mcell, false, true, 1e-5);
  }
}

TEST_CASE("non-collinear antiferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  auto const spins =
      noncollinear_tensors({Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, -1.0)});
  MagneticCell const mcell(cell, SiteTensors{spins});
  check(mcell, true, true, 1e-5);
}
