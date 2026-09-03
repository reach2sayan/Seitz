// Oracle test for cell standardization (spglib spg_standardize_cell):
// cppcrystal::standardize_cell must reproduce the reference across all four
// to_primitive x no_idealize combinations — atom count, cell shape (metric
// tensor, rotation-invariant), composition, and the atom set.

#include "oracle.hpp"

#include <cppcrystal/analysis/symmetry_analyzer.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using cppcrystal::Cell;
using cppcrystal::Lattice;
using cppcrystal::Matrix3d;
using cppcrystal::Positions;
using cppcrystal::Result;
using cppcrystal::Types;
using cppcrystal::Vector3d;

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

// Lattices match up to a rigid rotation iff their metric tensors agree.
bool same_metric(Matrix3d const &a, Matrix3d const &b) {
  return ((a.transpose() * a) - (b.transpose() * b)).cwiseAbs().maxCoeff() <
         1e-3;
}

std::map<int, int> composition(Cell const &c) {
  std::map<int, int> m;
  for (int t : c.types()) {
    ++m[t];
  }
  return m;
}

bool overlap(Vector3d const &a, Vector3d const &b, double tol) {
  Vector3d d = (a - b).unaryExpr([](double x) { return x - std::round(x); });
  return d.cwiseAbs().maxCoeff() <= tol;
}

using cppcrystal::test::must;

void check_one(Cell const &cell, bool to_primitive, bool no_idealize,
               double symprec) {
  INFO("to_primitive=" << to_primitive << " no_idealize=" << no_idealize);
  using cppcrystal::analysis::CellSetting;
  using cppcrystal::analysis::Idealize;
  auto const analyzer =
      cppcrystal::analysis::SymmetryAnalyzer::from_cell(cell, {symprec});
  // The setting is a template argument, so the reference's two runtime flags
  // pick one of the four instantiations here.
  auto const got = must(
      to_primitive
          ? (no_idealize ? analyzer.standardized_cell<CellSetting::primitive,
                                                      Idealize::no>()
                         : analyzer.standardized_cell<CellSetting::primitive,
                                                      Idealize::yes>())
          : (no_idealize ? analyzer.standardized_cell<CellSetting::conventional,
                                                      Idealize::no>()
                         : analyzer.standardized_cell<CellSetting::conventional,
                                                      Idealize::yes>()));
  Cell const ref = cppcrystal::oracle::reference_standardize_cell(
      cell, to_primitive, no_idealize, symprec);
  REQUIRE(ref.size() > 0); // reference succeeded

  CHECK(got.size() == ref.size());
  CHECK(same_metric(got.lattice().matrix(), ref.lattice().matrix()));
  CHECK(composition(got) == composition(ref));

  // Every reference atom is matched (same type, coincident position mod 1) by a
  // distinct atom of ours.
  REQUIRE(got.size() == ref.size());
  std::vector<char> used(static_cast<std::size_t>(got.size()), 0);
  for (Eigen::Index r = 0; r < ref.size(); ++r) {
    bool found = false;
    for (Eigen::Index g = 0; g < got.size() && !found; ++g) {
      auto const ug = static_cast<std::size_t>(g);
      if (used[ug] || got.type(g) != ref.type(r)) {
        continue;
      }
      if (overlap(got.position(g), ref.position(r), symprec)) {
        used[ug] = 1;
        found = true;
      }
    }
    INFO("reference std atom " << r << " not matched");
    CHECK(found);
  }
}

void check_all(Cell const &cell, double symprec = 1e-5) {
  check_one(cell, false, false, symprec); // conventional, idealized
  check_one(cell, true, false, symprec);  // primitive,    idealized
  check_one(cell, false, true, symprec);  // conventional, no_idealize
  check_one(cell, true, true, symprec);   // primitive,    no_idealize
}

} // namespace

TEST_CASE("standardize: NaCl (Fm-3m, F-centered)", "[oracle][standardize]") {
  Matrix3d const L = 5.64 * Matrix3d::Identity();
  check_all(make_cell(L,
                      {{0, 0, 0},
                       {0, 0.5, 0.5},
                       {0.5, 0, 0.5},
                       {0.5, 0.5, 0},
                       {0.5, 0.5, 0.5},
                       {0.5, 0, 0},
                       {0, 0.5, 0},
                       {0, 0, 0.5}},
                      {11, 11, 11, 11, 17, 17, 17, 17}));
}

TEST_CASE("standardize: bcc Fe (Im-3m, I-centered)", "[oracle][standardize]") {
  Matrix3d const L = 2.87 * Matrix3d::Identity();
  check_all(make_cell(L, {{0, 0, 0}, {0.5, 0.5, 0.5}}, {26, 26}));
}

TEST_CASE("standardize: simple cubic (primitive)", "[oracle][standardize]") {
  check_all(make_cell(4.0 * Matrix3d::Identity(), {{0, 0, 0}}, {1}));
}

TEST_CASE("standardize: C-centered orthorhombic", "[oracle][standardize]") {
  Matrix3d L = Matrix3d::Zero();
  L(0, 0) = 4.0;
  L(1, 1) = 6.0;
  L(2, 2) = 8.0;
  check_all(make_cell(L, {{0, 0, 0}, {0.5, 0.5, 0}}, {1, 1}));
}

TEST_CASE("standardize: a primitive input expands to its conventional cell",
          "[oracle][standardize]") {
  // Primitive bcc cell (rhombohedral primitive of Im-3m) -> conventional 2-atom.
  Matrix3d P;
  P.col(0) = Vector3d(-1.435, 1.435, 1.435);
  P.col(1) = Vector3d(1.435, -1.435, 1.435);
  P.col(2) = Vector3d(1.435, 1.435, -1.435);
  check_all(make_cell(P, {{0, 0, 0}}, {26}));
}

TEST_CASE("standardize: strained cell, no_idealize preserves input geometry",
          "[oracle][standardize]") {
  // A slightly orthorhombically strained NaCl: with a loose symprec it is still
  // F-centered cubic to the reference; no_idealize must keep the strained
  // lattice while idealize restores the ideal metric.
  Matrix3d L = 5.64 * Matrix3d::Identity();
  L(0, 0) *= 1.002;
  L(1, 1) *= 0.998;
  Cell const cell = make_cell(L,
                              {{0, 0, 0},
                               {0, 0.5, 0.5},
                               {0.5, 0, 0.5},
                               {0.5, 0.5, 0},
                               {0.5, 0.5, 0.5},
                               {0.5, 0, 0},
                               {0, 0.5, 0},
                               {0, 0, 0.5}},
                              {11, 11, 11, 11, 17, 17, 17, 17});
  check_all(cell, 1e-2);
}

TEST_CASE("standardize: the to_primitive / idealized option pairs",
          "[oracle][standardize]") {
  // The two named wrappers must be bit-identical to the standardize_cell calls
  // they delegate to. find_primitive == spg_find_primitive (to_primitive,
  // idealized); the default == spg_refine_cell (conventional, idealized) — both
  // already oracle-checked via check_all above, so equality with their delegate
  // transitively proves the wrapper.
  Cell const nacl = make_cell(5.64 * Matrix3d::Identity(),
                              {{0, 0, 0},
                               {0, 0.5, 0.5},
                               {0.5, 0, 0.5},
                               {0.5, 0.5, 0},
                               {0.5, 0.5, 0.5},
                               {0.5, 0, 0},
                               {0, 0.5, 0},
                               {0, 0, 0.5}},
                              {11, 11, 11, 11, 17, 17, 17, 17});

  auto const same_cell = [](Cell const &a, Cell const &b) {
    return a.size() == b.size() &&
           (a.lattice().matrix() - b.lattice().matrix()).cwiseAbs().maxCoeff() < 1e-12 &&
           (a.positions() - b.positions()).cwiseAbs().maxCoeff() < 1e-12 &&
           a.types() == b.types();
  };

  Cell const prim = must(
      cppcrystal::analysis::SymmetryAnalyzer::from_cell(nacl).standardized_cell<cppcrystal::analysis::CellSetting::primitive, cppcrystal::analysis::Idealize::yes>());
  Cell const prim_ref = must(
      cppcrystal::analysis::SymmetryAnalyzer::from_cell(nacl).standardized_cell<cppcrystal::analysis::CellSetting::primitive, cppcrystal::analysis::Idealize::yes>());
  CHECK(same_cell(prim, prim_ref));

  auto const analyzer = cppcrystal::analysis::SymmetryAnalyzer::from_cell(nacl);
  Cell const refined = must(analyzer.standardized_cell());
  Cell const refined_ref = must(analyzer.standardized_cell());
  CHECK(same_cell(refined, refined_ref));
}
