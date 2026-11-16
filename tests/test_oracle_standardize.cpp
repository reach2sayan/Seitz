// Oracle test for cell standardization (spglib spg_standardize_cell):
// seitz::standardize_cell must reproduce the reference across all four
// to_primitive x no_idealize combinations — atom count, cell shape (metric
// tensor, rotation-invariant), composition, and the atom set.

#include "oracle.hpp"
#include "structures.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>

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

using seitz::test::make_cell;

using seitz::Cell;
using seitz::Lattice;
using seitz::Matrix3d;
using seitz::Positions;
using seitz::Result;
using seitz::Types;
using seitz::Vector3d;

std::map<int, int> composition(Cell const &c) {
  std::map<int, int> m;
  for (int t : c.types()) {
    ++m[t];
  }
  return m;
}

using seitz::test::must;

void check_one(Cell const &cell, bool to_primitive, bool no_idealize,
               double symprec) {
  INFO("to_primitive=" << to_primitive << " no_idealize=" << no_idealize);
  using seitz::analysis::CellSetting;
  using seitz::analysis::Idealize;
  auto const analyzer =
      seitz::analysis::SymmetryAnalyzer::from_cell(cell, {symprec});
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
  Cell const ref = seitz::oracle::reference_standardize_cell(
      cell, to_primitive, no_idealize, symprec);
  REQUIRE(ref.size() > 0); // reference succeeded

  CHECK(got.size() == ref.size());
  CHECK(seitz::test::same_metric(got.lattice().matrix(), ref.lattice().matrix(),
                                 1e-3));
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
      if (seitz::test::coincide_mod1(got.position(g), ref.position(r),
                                     symprec)) {
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
  check_all(seitz::test::nacl());
}

TEST_CASE("standardize: bcc Fe (Im-3m, I-centered)", "[oracle][standardize]") {
  Matrix3d const L = 2.87 * Matrix3d::Identity();
  check_all(make_cell(L, {{0, 0, 0}, {0.5, 0.5, 0.5}}, {26, 26}));
}

TEST_CASE("standardize: simple cubic (primitive)", "[oracle][standardize]") {
  check_all(make_cell(4.0 * Matrix3d::Identity(), {{0, 0, 0}}, {1}));
}

TEST_CASE("standardize: C-centered orthorhombic", "[oracle][standardize]") {
  check_all(make_cell(seitz::test::orthorhombic(4.0, 6.0, 8.0),
                      {{0, 0, 0}, {0.5, 0.5, 0}}, {1, 1}));
}

TEST_CASE("standardize: a primitive input expands to its conventional cell",
          "[oracle][standardize]") {
  // Primitive bcc cell (rhombohedral primitive of Im-3m) -> conventional
  // 2-atom.
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
  Cell const cell = seitz::test::nacl().with_lattice(Lattice{L});
  check_all(cell, 1e-2);
}
