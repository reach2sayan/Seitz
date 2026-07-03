// Non-oracle tests for standardize_cell: round-trip checks against the
// determination pipeline, without the reference spglib.
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/standardize.hpp>

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <utility>

using namespace cppcrystal;

namespace {
template <class T> T must(Result<T> r) {
  return leaf::try_handle_all(
      [&]() -> Result<T> { return std::move(r); },
      [](leaf::error_info const &) -> T {
        FAIL("unexpected error");
        throw std::logic_error("unreachable");
      });
}

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
  return Cell(L, p, Types{11, 11, 11, 11, 17, 17, 17, 17});
}

int spacegroup_of(Cell const &c) {
  return must(analysis::SymmetryAnalyzer::from_cell(c).spacegroup_number());
}
} // namespace

TEST_CASE("standardize_cell: conventional vs primitive atom counts",
          "[standardize]") {
  Cell const nacl = nacl_conventional();
  // F-centering: conventional has 4x the primitive atoms.
  CHECK(must(standardize_cell(nacl, {.to_primitive = false})).size() == 8);
  CHECK(must(standardize_cell(nacl, {.to_primitive = true})).size() == 2);
  // All settings describe the same space group.
  CHECK(spacegroup_of(must(standardize_cell(nacl, {false, false}))) == 225);
  CHECK(spacegroup_of(must(standardize_cell(nacl, {true, false}))) == 225);
  CHECK(spacegroup_of(must(standardize_cell(nacl, {false, true}))) == 225);
  CHECK(spacegroup_of(must(standardize_cell(nacl, {true, true}))) == 225);
}

TEST_CASE("standardize_cell: no_idealize keeps an undistorted cell unchanged",
          "[standardize]") {
  // For an already-ideal cell, the idealized and non-idealized conventional
  // standardizations agree on the lattice metric.
  Cell const nacl = nacl_conventional();
  auto const ideal = must(standardize_cell(nacl, {false, false}));
  auto const raw = must(standardize_cell(nacl, {false, true}));
  Matrix3d const gi = ideal.lattice().transpose() * ideal.lattice();
  Matrix3d const gr = raw.lattice().transpose() * raw.lattice();
  CHECK((gi - gr).cwiseAbs().maxCoeff() < 1e-6);
}

TEST_CASE("standardize_cell: no_idealize preserves a strained input lattice",
          "[standardize]") {
  // Strain the cell so its true metric is non-cubic but, within a loose
  // symprec (0.1 A > the ~0.056 A length spread), still Fm-3m. no_idealize must
  // keep the strained vectors; idealize must restore a == b == c.
  Matrix3d L = 5.64 * Matrix3d::Identity();
  L(0, 0) *= 1.005;
  L(2, 2) *= 0.995;
  Cell strained = nacl_conventional();
  strained.lattice() = L;

  auto const raw = must(standardize_cell(strained, {false, true}, 1e-1));
  // The conventional no_idealize lattice keeps the input's distinct a/b/c
  // lengths (volume and column norms preserved up to centering reordering),
  // unlike the idealized cell which would force a == b == c.
  double const va = raw.lattice().col(0).norm();
  double const vb = raw.lattice().col(1).norm();
  double const vc = raw.lattice().col(2).norm();
  CHECK(std::abs(va - vc) > 1e-3); // genuinely non-cubic, distortion retained

  auto const ideal = must(standardize_cell(strained, {false, false}, 1e-1));
  CHECK(std::abs(ideal.lattice().col(0).norm() -
                 ideal.lattice().col(2).norm()) < 1e-3); // idealized to cubic
}

TEST_CASE("refine_cell equals the default standardize_cell", "[standardize]") {
  Cell const nacl = nacl_conventional();
  auto const a = must(refine_cell(nacl));
  auto const b = must(standardize_cell(nacl, {}));
  REQUIRE(a.size() == b.size());
  CHECK(a.lattice().isApprox(b.lattice()));
}
