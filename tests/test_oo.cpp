// Tests for the object-oriented layer: the analysis facades (SymmetryAnalyzer),
// the standalone group objects (SpaceGroup / WyckoffPosition), and crystal
// generation. The group tests assert the orbit-stabilizer invariant across all
// 230 space groups; the generation test round-trips a generated cell back
// through the analyzer.
#include <spglib/analysis/symmetry_analyzer.hpp>
#include <spglib/core/overlap.hpp>
#include <spglib/dataset.hpp>
#include <spglib/generate/crystal_builder.hpp>
#include <spglib/generate/distance_check.hpp>
#include <spglib/group/space_group.hpp>
#include <spglib/group/subgroup_graph.hpp>

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

  auto gen = must(generate::random_crystal(sg, comp, {.seed = 42u}));
  REQUIRE(gen.cell.size() == 8);
  std::set<int> kinds(gen.cell.types().begin(), gen.cell.types().end());
  REQUIRE(kinds == std::set<int>{11, 17});

  // The assembled structure actually has the requested symmetry.
  auto analyzer = analysis::SymmetryAnalyzer::from_cell(gen.cell);
  REQUIRE(must(analyzer.spacegroup_number()) == 225);
}

TEST_CASE("generated structures recover their target space group", "[generate]") {
  // A spread of targets with fixed-position orbits, each round-tripped through
  // the determination pipeline (self-validation, no external oracle).
  struct Case {
    int number;
    generate::Composition comp;
  };
  auto const cases = {
      Case{221, {{1, 1}}},          // Pm-3m, 1a
      Case{229, {{26, 2}}},         // Im-3m (bcc Fe), 2a
      Case{225, {{11, 4}, {17, 4}}} // Fm-3m (NaCl)
  };
  for (auto const &c : cases) {
    auto sg = must(group::SpaceGroup::from_number(c.number));
    auto gen = must(generate::random_crystal(sg, c.comp, {.seed = 7u}));
    auto analyzer = analysis::SymmetryAnalyzer::from_cell(gen.cell);
    REQUIRE(must(analyzer.spacegroup_number()) == c.number);
  }
}

TEST_CASE("generated structures are free of interatomic clashes", "[generate]") {
  auto sg = must(group::SpaceGroup::from_number(225));
  generate::Composition const comp{{11, 4}, {17, 4}};
  auto gen = must(generate::random_crystal(sg, comp, {.seed = 99u}));
  // The builder only returns distance-valid cells; confirm independently.
  REQUIRE(generate::distances_valid(gen.cell));
}

TEST_CASE("generation is deterministic in the seed", "[generate]") {
  auto sg = must(group::SpaceGroup::from_number(225));
  generate::Composition const comp{{11, 4}, {17, 4}};
  auto a = must(generate::random_crystal(sg, comp, {.seed = 123u}));
  auto b = must(generate::random_crystal(sg, comp, {.seed = 123u}));
  REQUIRE(a.cell.positions().isApprox(b.cell.positions()));
  REQUIRE(a.cell.lattice().isApprox(b.cell.lattice()));
}

TEST_CASE("incompatible composition is rejected", "[generate]") {
  auto sg = must(group::SpaceGroup::from_number(225)); // Fm-3m, smallest mult 4
  // Three atoms cannot fill any combination of multiplicity-4-or-more orbits.
  REQUIRE_FALSE(generate::check_compatible(sg, generate::Composition{{1, 3}}));
}

TEST_CASE("orbit-stabilizer invariant holds for all 80 layer groups",
          "[layergen]") {
  for (int number = 1; number <= 80; ++number) {
    auto lg = must(group::SpaceGroup::from_layer_number(number));
    REQUIRE(lg.number() == number);
    REQUIRE(lg.hall_number() < 0); // negative-Hall (layer) convention
    auto const nops = static_cast<int>(lg.operations().size());
    REQUIRE(nops > 0);
    for (auto const &wp : lg.wyckoffs()) {
      REQUIRE(wp.multiplicity() * static_cast<int>(wp.operations().size()) ==
              nops);
    }
  }
}

TEST_CASE("generated layer structures carry their full layer symmetry",
          "[layergen]") {
  // The direct self-validation: a generated 2D-periodic structure must be
  // invariant under EVERY operation of its layer group (aperiodic-aware overlap).
  // This covers all crystal systems and the c-flipping groups (inversion,
  // in-plane 2-folds/mirrors, horizontal mirrors), which the orbit must build
  // without folding the non-periodic axis.
  for (int number : {1, 2, 19, 49, 61, 65}) {
    INFO("layer group " << number);
    auto lg = must(group::SpaceGroup::from_layer_number(number));
    int const m = lg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::random_layer_crystal(
        lg, comp,
        {.volume_factor = 4.0, .seed = 13u, .general_position_only = true}));
    REQUIRE(gen.cell.aperiodic_axis() == 2);
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(generate::distances_valid(gen.cell));

    OverlapChecker checker(gen.cell);
    for (auto const &op : lg.operations()) {
      REQUIRE(checker.check_total_overlap(op.translation, op.rotation, 1e-3));
    }
  }
}

TEST_CASE("layer crystal round-trips through get_layer_dataset", "[layergen]") {
  // End-to-end self-validation through the layer determination. Restricted to
  // c-preserving layer groups (one atomic plane): get_layer_dataset currently
  // does not resolve multi-level (c-flipping) layers, so those are covered by
  // the direct op-invariance test above rather than a determination round-trip.
  for (int number : {1, 49, 55, 65}) {
    INFO("layer group " << number);
    auto lg = must(group::SpaceGroup::from_layer_number(number));
    int const m = lg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::random_layer_crystal(
        lg, comp,
        {.volume_factor = 4.0, .seed = 13u, .general_position_only = true}));
    auto ds = must(get_layer_dataset(gen.cell, 2, 1e-4));
    REQUIRE(ds.spacegroup_number == number); // exact layer-group recovery
  }
}

namespace {
// Point-group order of a space group: distinct rotation parts of its operations.
int pg_order(int number) {
  auto sg = must(group::SpaceGroup::from_number(number));
  std::vector<Matrix3i> rots;
  for (auto const &op : sg.operations()) {
    if (std::ranges::none_of(
            rots, [&](Matrix3i const &r) { return r == op.rotation; })) {
      rots.push_back(op.rotation);
    }
  }
  return static_cast<int>(rots.size());
}

bool has_relation(std::vector<group::SubgroupRelation> const &rels, int number,
                  int index) {
  return std::ranges::any_of(rels, [&](group::SubgroupRelation const &r) {
    return r.number == number && r.index == index;
  });
}
} // namespace

TEST_CASE("Pm-3m maximal t-subgroups match the textbook relations",
          "[subgroup]") {
  auto const &graph = group::SubgroupGraph::instance();
  auto const subs = graph.maximal_subgroups(221); // Oh
  // The five maximal translationengleiche subgroups of Pm-3m (ITA): Pm-3 (200,
  // Th), P432 (207, O), P-43m (215, Td) at index 2; P4/mmm (123, D4h) at index
  // 3; R-3m (166, D3d, along a body diagonal) at index 4 — maximal subgroups of
  // a non-nilpotent group may have composite index.
  REQUIRE(has_relation(subs, 200, 2));
  REQUIRE(has_relation(subs, 207, 2));
  REQUIRE(has_relation(subs, 215, 2));
  REQUIRE(has_relation(subs, 123, 3));
  REQUIRE(has_relation(subs, 166, 4));
  REQUIRE(subs.size() == 5);
}

TEST_CASE("P-1 has the single maximal t-subgroup P1", "[subgroup]") {
  auto const &graph = group::SubgroupGraph::instance();
  auto const subs = graph.maximal_subgroups(2);
  REQUIRE(subs.size() == 1);
  REQUIRE(subs.front().number == 1);
  REQUIRE(subs.front().index == 2);
  // P1 itself has no proper subgroup.
  REQUIRE(graph.maximal_subgroups(1).empty());
}

TEST_CASE("every t-subgroup edge is order-consistent", "[subgroup]") {
  auto const &graph = group::SubgroupGraph::instance();
  int edges = 0;
  for (int n = 1; n <= 230; ++n) {
    for (auto const &rel : graph.maximal_subgroups(n)) {
      REQUIRE(rel.index >= 2);
      REQUIRE(rel.number != n);
      // A translationengleiche subgroup keeps the lattice, so the index equals
      // the ratio of point-group orders.
      REQUIRE(pg_order(n) == rel.index * pg_order(rel.number));
      ++edges;
    }
  }
  REQUIRE(edges > 200); // the graph is richly connected
}

TEST_CASE("reachability and symmetry-breaking paths", "[subgroup]") {
  auto const &graph = group::SubgroupGraph::instance();
  REQUIRE(graph.is_subgroup(221, 221)); // reflexive
  REQUIRE(graph.is_subgroup(1, 221));   // P1 is a t-subgroup of Pm-3m
  REQUIRE_FALSE(graph.is_subgroup(221, 1));

  auto chain = graph.path(221, 1);
  REQUIRE(chain.has_value());
  REQUIRE(chain->front() == 221);
  REQUIRE(chain->back() == 1);
  // Each step descends along a real maximal-subgroup edge (any index).
  for (std::size_t i = 0; i + 1 < chain->size(); ++i) {
    auto const subs = graph.maximal_subgroups((*chain)[i]);
    REQUIRE(std::ranges::any_of(subs, [&](group::SubgroupRelation const &r) {
      return r.number == (*chain)[i + 1];
    }));
  }
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
