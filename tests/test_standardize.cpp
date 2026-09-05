// Non-oracle tests for standardize_cell: round-trip checks against the
// determination pipeline, without the reference spglib.
#include <seitz/analysis/symmetry_analyzer.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <utility>

using namespace seitz;

namespace {
using seitz::test::must;

Cell nacl_conventional() {
  Matrix3d const L = 5.64 * Matrix3d::Identity();
  Positions p(8, 3);
  p.row(0) << 0, 0, 0;
  p.row(1) << 0, 0.5, 0.5;
  p.row(2) << 0.5, 0, 0.5;
  p.row(3) << 0.5, 0.5, 0;
  p.row(4) << 0.5, 0.5, 0.5;
  p.row(5) << 0.5, 0, 0;
  p.row(6) << 0, 0.5, 0;
  p.row(7) << 0, 0, 0.5;
  return Cell(Lattice{L}, p, Types{11, 11, 11, 11, 17, 17, 17, 17});
}

int spacegroup_of(Cell const &c) {
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(c);
  return data::spacegroup_type(must(analyzer.hall())).number;
}
} // namespace

TEST_CASE("standardize_cell: conventional vs primitive atom counts",
          "[standardize]") {
  Cell const nacl = nacl_conventional();
  // F-centering: conventional has 4x the primitive atoms.
  CHECK(must(analysis::SymmetryAnalyzer::from_cell(nacl)
                 .standardized_cell<analysis::CellSetting::conventional,
                                    analysis::Idealize::yes>())
            .size() == 8);
  CHECK(must(analysis::SymmetryAnalyzer::from_cell(nacl)
                 .standardized_cell<analysis::CellSetting::primitive,
                                    analysis::Idealize::yes>())
            .size() == 2);
  // All four settings describe the same space group. The setting is a template
  // argument now, so the four combinations are named rather than looped.
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(nacl);
  using analysis::CellSetting;
  using analysis::Idealize;
  CHECK(spacegroup_of(must(analyzer.standardized_cell<CellSetting::conventional,
                                                      Idealize::yes>())) ==
        225);
  CHECK(spacegroup_of(must(analyzer.standardized_cell<CellSetting::conventional,
                                                      Idealize::no>())) == 225);
  CHECK(spacegroup_of(must(analyzer.standardized_cell<CellSetting::primitive,
                                                      Idealize::yes>())) ==
        225);
  CHECK(spacegroup_of(must(analyzer.standardized_cell<CellSetting::primitive,
                                                      Idealize::no>())) == 225);
}

TEST_CASE("standardize_cell: Idealize::no keeps an undistorted cell unchanged",
          "[standardize]") {
  // For an already-ideal cell, the idealized and non-idealized conventional
  // standardizations agree on the lattice metric.
  Cell const nacl = nacl_conventional();
  auto const ideal =
      must(analysis::SymmetryAnalyzer::from_cell(nacl)
               .standardized_cell<analysis::CellSetting::conventional,
                                  analysis::Idealize::yes>());
  auto const raw =
      must(analysis::SymmetryAnalyzer::from_cell(nacl)
               .standardized_cell<analysis::CellSetting::conventional,
                                  analysis::Idealize::no>());
  Matrix3d const gi =
      ideal.lattice().matrix().transpose() * ideal.lattice().matrix();
  Matrix3d const gr =
      raw.lattice().matrix().transpose() * raw.lattice().matrix();
  CHECK((gi - gr).cwiseAbs().maxCoeff() < 1e-6);
}

TEST_CASE("standardize_cell: no_idealize preserves a strained input lattice",
          "[standardize]") {
  // Strain the cell so its true metric is non-cubic but, within a loose
  // symprec (0.1 A > the ~0.056 A length spread), still Fm-3m. Idealize::no
  // must keep the strained vectors; Idealize::yes must restore a == b == c.
  Matrix3d L = 5.64 * Matrix3d::Identity();
  L(0, 0) *= 1.005;
  L(2, 2) *= 0.995;
  Cell const strained = nacl_conventional().with_lattice(Lattice{L});

  auto const raw =
      must(analysis::SymmetryAnalyzer::from_cell(strained, {1e-1})
               .standardized_cell<analysis::CellSetting::conventional,
                                  analysis::Idealize::no>());
  // The conventional Idealize::no lattice keeps the input's distinct a/b/c
  // lengths (volume and column norms preserved up to centering reordering),
  // unlike the idealized cell which would force a == b == c.
  double const va = raw.lattice().matrix().col(0).norm();
  double const vb = raw.lattice().matrix().col(1).norm();
  double const vc = raw.lattice().matrix().col(2).norm();
  CHECK(std::abs(va - vc) > 1e-3); // genuinely non-cubic, distortion retained

  auto const ideal =
      must(analysis::SymmetryAnalyzer::from_cell(strained, {1e-1})
               .standardized_cell<analysis::CellSetting::conventional,
                                  analysis::Idealize::yes>());
  CHECK(std::abs(ideal.lattice().matrix().col(0).norm() -
                 ideal.lattice().matrix().col(2).norm()) <
        1e-3); // idealized to cubic
}
