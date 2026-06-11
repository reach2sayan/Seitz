// Tests for the object-oriented layer: the analysis facades (SymmetryAnalyzer),
// the standalone group objects (SpaceGroup / WyckoffPosition), and crystal
// generation. The group tests assert the orbit-stabilizer invariant across all
// 230 space groups; the generation test round-trips a generated cell back
// through the analyzer.
#include <spglib/analysis/symmetry_analyzer.hpp>
#include <spglib/generate/crystal_builder.hpp>
#include <spglib/group/space_group.hpp>

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <set>
#include <stdexcept>
#include <utility>

using namespace spglib;

namespace {
// Unwrap a move-only Result<T>, failing the test on an error result.
template <class T> T must(Result<T> r) {
  return leaf::try_handle_all(
      [&]() -> Result<T> { return std::move(r); },
      [](leaf::error_info const &) -> T {
        FAIL("unexpected error result");
        throw std::logic_error("unreachable");
      });
}
} // namespace

TEST_CASE("orbit-stabilizer invariant holds for all 230 space groups",
          "[group]") {
  int checked = 0;
  for (int number = 1; number <= 230; ++number) {
    auto sg = must(group::SpaceGroup::from_number(number));
    auto const nops = static_cast<int>(sg.operations().size());
    REQUIRE(nops > 0);
    for (auto const &wp : sg.wyckoffs()) {
      auto const order = static_cast<int>(wp.operations().size());
      // multiplicity * |site-symmetry group| == |conventional operations|
      REQUIRE(wp.multiplicity() * order == nops);
      ++checked;
    }
    // The general position is last, fully free, with trivial site symmetry.
    auto const &general = sg.wyckoffs().back();
    REQUIRE(general.degrees_of_freedom() == 3);
    REQUIRE(general.operations().size() == 1);
  }
  REQUIRE(checked > 1000);
}

TEST_CASE("SpaceGroup Pm-3m (221) matches ITA reference", "[group]") {
  auto sg = must(group::SpaceGroup::from_number(221));
  REQUIRE(sg.number() == 221);
  REQUIRE(sg.operations().size() == 48);

  auto a = must(sg.wyckoff('a'));
  REQUIRE(a->multiplicity() == 1);
  REQUIRE(a->degrees_of_freedom() == 0);
  REQUIRE(a->site_symmetry() == "m-3m");

  // 1a has no free coordinate: its orbit is the single fixed point.
  Positions const orbit = a->get_all_positions(Vector3d{0.3, 0.4, 0.5});
  REQUIRE(orbit.rows() == 1);

  REQUIRE(sg.wyckoffs().back().multiplicity() == 48);
}

TEST_CASE("WyckoffPosition orbit expansion respects multiplicity", "[group]") {
  auto sg = must(group::SpaceGroup::from_number(225)); // Fm-3m
  auto const &general = sg.wyckoffs().back();
  Positions const orbit = general.get_all_positions(Vector3d{0.11, 0.23, 0.37});
  REQUIRE(orbit.rows() == general.multiplicity());
}

TEST_CASE("from_number rejects out-of-range numbers", "[group]") {
  bool errored = leaf::try_handle_all(
      [&]() -> Result<bool> {
        BOOST_LEAF_AUTO(sg, group::SpaceGroup::from_number(231));
        (void)sg;
        return false;
      },
      [](leaf::error_info const &) { return true; });
  REQUIRE(errored);
}

TEST_CASE("crystal generation round-trips through the analyzer", "[generate]") {
  auto sg = must(group::SpaceGroup::from_number(225)); // Fm-3m
  // NaCl: Na (type 11) and Cl (type 17), four of each (the 4a / 4b orbits).
  generate::Composition const comp{{11, 4}, {17, 4}};
  REQUIRE(generate::check_compatible(sg, comp));

  auto cell = must(generate::random_crystal(sg, comp, 1.0, 42u));
  REQUIRE(cell.size() == 8);
  std::set<int> kinds(cell.types().begin(), cell.types().end());
  REQUIRE(kinds == std::set<int>{11, 17});

  auto analyzer = analysis::SymmetryAnalyzer::from_cell(cell);
  REQUIRE(must(analyzer.spacegroup_number()) == 225);
}

TEST_CASE("incompatible composition is rejected", "[generate]") {
  auto sg = must(group::SpaceGroup::from_number(225)); // Fm-3m, smallest mult 4
  // Three atoms cannot fill any combination of multiplicity-4-or-more orbits.
  REQUIRE_FALSE(generate::check_compatible(sg, generate::Composition{{1, 3}}));
}

TEST_CASE("SymmetryAnalyzer memoizes a consistent dataset", "[analysis]") {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  Cell const cell{Matrix3d::Identity(), pos, Types{1}};

  auto analyzer = analysis::SymmetryAnalyzer::from_cell(cell);
  int const first = must(analyzer.spacegroup_number());
  int const cached = must(analyzer.spacegroup_number());
  REQUIRE(first == 221); // simple cubic, one atom -> Pm-3m
  REQUIRE(first == cached);

  // Projections agree with the dataset.
  auto ds = must(analyzer.dataset());
  REQUIRE(ds.spacegroup_number == first);
  REQUIRE(must(analyzer.operations()).size() == ds.operations.size());
}
