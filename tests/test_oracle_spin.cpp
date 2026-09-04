// Oracle test for the magnetic symmetry search (Phase 7, spin.c): the magnetic
// operations / equivalent atoms / primitive lattice found by
// spin::operations_with_site_tensors must match the reference
// spg_get_symmetry_with_site_tensors.

#include "oracle.hpp"

#include "spin/search.hpp"
#include "symmetry/search.hpp"
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/operation_set.hpp>

#include "math/integer_matrix.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::CollinearTensors;
using cppcrystal::Lattice;
using cppcrystal::MagneticCell;
using cppcrystal::MagneticOperations;
using cppcrystal::Matrix3d;
using cppcrystal::Matrix3i;
using cppcrystal::noncollinear_tensors;
using cppcrystal::Positions;
using cppcrystal::SiteTensors;
using cppcrystal::Vector3d;

Cell make_cell(double a, std::vector<std::array<double, 3>> const &pos,
               std::vector<int> const &types) {
  Matrix3d lattice = a * Matrix3d::Identity();
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(Lattice{lattice}, p, types);
}

// Reference magnetic symmetry via spg_get_symmetry_with_site_tensors. `tensors`
// is flat: collinear (rank 0) = one per atom, non-collinear (rank 1) = three
// per atom. Returns operations with time_reversal derived from spin_flips
// (spin_flips = 1 - 2*timerev, so timerev <=> spin_flips == -1).
struct Reference {
  std::vector<cppcrystal::MagneticSymmetryOperation> operations;
  std::vector<int> equivalent_atoms;
  Matrix3d primitive_lattice;
};

Reference reference_magnetic(Cell const &cell,
                             std::vector<double> const &tensors,
                             int tensor_rank, bool with_time_reversal,
                             bool is_axial, double symprec) {
  cppcrystal::oracle::CCell c(cell);
  int const n = c.num_atom();
  int const max_size = 384;
  std::vector<int> rot(static_cast<std::size_t>(9 * max_size));
  std::vector<double> trans(static_cast<std::size_t>(3 * max_size));
  std::vector<int> equiv(static_cast<std::size_t>(n));
  std::vector<int> spin_flips(static_cast<std::size_t>(max_size));
  double prim[3][3];

  int const num = spg_get_symmetry_with_site_tensors(
      reinterpret_cast<int(*)[3][3]>(rot.data()),
      reinterpret_cast<double(*)[3]>(trans.data()), equiv.data(), prim,
      spin_flips.data(), max_size, c.lattice, c.pos(), c.types.data(),
      tensors.data(), tensor_rank, n, with_time_reversal, is_axial, symprec);
  REQUIRE(num > 0);

  Reference ref;
  ref.operations.reserve(static_cast<std::size_t>(num));
  for (int s = 0; s < num; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = rot[static_cast<std::size_t>(9 * s + 3 * i + j)];
    Vector3d t(trans[static_cast<std::size_t>(3 * s + 0)],
               trans[static_cast<std::size_t>(3 * s + 1)],
               trans[static_cast<std::size_t>(3 * s + 2)]);
    bool const time_reversal = spin_flips[static_cast<std::size_t>(s)] == -1;
    ref.operations.push_back({{r, t}, time_reversal});
  }
  ref.equivalent_atoms = equiv;
  cppcrystal::oracle::from_c_lattice(ref.primitive_lattice, prim);
  return ref;
}

bool contains_operation(auto const &ops,
                        cppcrystal::MagneticSymmetryOperation const &target,
                        double symprec) {
  for (auto const &op : ops) {
    if (op.spatial.rotation != target.spatial.rotation ||
        op.time_reversal != target.time_reversal) {
      continue;
    }
    Vector3d d = op.spatial.translation - target.spatial.translation;
    d = d.unaryExpr([](double x) { return x - std::round(x); });
    if (d.cwiseAbs().maxCoeff() <= symprec) {
      return true;
    }
  }
  return false;
}

// Both lattices describe the same lattice iff L_a^-1 . L_b is an integer
// unimodular matrix (a change of primitive basis).
bool same_lattice(Matrix3d const &a, Matrix3d const &b) {
  Matrix3d const rel = a.inverse() * b;
  Matrix3i const rounded = cppcrystal::math::round_to_int(rel);
  if ((rel - rounded.cast<double>()).cwiseAbs().maxCoeff() > 1e-6) {
    return false;
  }
  return std::abs(rounded.determinant()) == 1;
}

void check(MagneticCell const &input, std::vector<double> const &tensors,
           int tensor_rank, bool with_time_reversal, bool is_axial,
           double symprec) {
  MagneticCell const mcell(input.cell(), input.tensors(),
                           is_axial ? cppcrystal::TensorKind::axial
                                    : cppcrystal::TensorKind::polar);
  cppcrystal::symmetry::SymmetrySearch<cppcrystal::GroupFamily::space> const
      search(mcell.cell(), {symprec});
  auto const sym_nonspin = search.operations();
  REQUIRE(sym_nonspin);

  cppcrystal::spin::SpinSearch const spin_search(mcell, sym_nonspin.value(),
                                                 {{symprec}});
  auto const got =
      with_time_reversal
          ? spin_search.operations<cppcrystal::TimeReversal::on>()
          : spin_search.operations<cppcrystal::TimeReversal::off>();
  REQUIRE(got);

  auto const ref = reference_magnetic(mcell.cell(), tensors, tensor_rank,
                                      with_time_reversal, is_axial, symprec);

  // Operations match as a set (the spatial-symmetry order may differ).
  REQUIRE(got->operations.size() == ref.operations.size());
  for (auto const &op : ref.operations) {
    INFO("missing reference operation");
    REQUIRE(contains_operation(got->operations, op, symprec));
  }
  REQUIRE(got->equivalent_atoms == ref.equivalent_atoms);
  REQUIRE(same_lattice(got->primitive_lattice, ref.primitive_lattice));
}

} // namespace

TEST_CASE("collinear ferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}}, {0});
  CollinearTensors const spins{1.0};
  MagneticCell const mcell(cell, SiteTensors{spins});
  check(mcell, {1.0}, 0, true, true, 1e-5);
}

TEST_CASE("collinear antiferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  CollinearTensors const spins{1.0, -1.0};
  MagneticCell const mcell(cell, SiteTensors{spins});
  SECTION("with time reversal (family group)") {
    check(mcell, {1.0, -1.0}, 0, true, true, 1e-5);
  }
  SECTION("without time reversal (maximal subgroup)") {
    check(mcell, {1.0, -1.0}, 0, false, true, 1e-5);
  }
}

TEST_CASE("non-collinear antiferromagnet magnetic symmetry", "[oracle][spin]") {
  Cell const cell = make_cell(4.0, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, {0, 0});
  auto const spins =
      noncollinear_tensors({Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, -1.0)});
  MagneticCell const mcell(cell, SiteTensors{spins});
  check(mcell, {0.0, 0.0, 1.0, 0.0, 0.0, -1.0}, 1, true, true, 1e-5);
}
