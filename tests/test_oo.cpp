// The object-oriented layer: the analysis facades, the standalone group
// objects, and generation. The group tests assert orbit-stabilizer across all
// 230 space groups; generation round-trips a cell back through the analyzer.
#include "core/overlap.hpp"
#include "data/rod_database.hpp"
#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/generate/distance_check.hpp>
#include <seitz/generate/generator.hpp>
#include <seitz/group/point_group.hpp>
#include <seitz/group/rod_group.hpp>
#include <seitz/group/space_group.hpp>
#include <seitz/warmup.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <ranges>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

using namespace seitz;

namespace {
// Unwrap a move-only Result<T>, failing the test on an error result.
using seitz::test::errored;
using seitz::test::must;

// Order of the point group a tabulated site-symmetry symbol names. These are
// oriented Hermann-Mauguin symbols ("..2", "m.mm"): dots mark direction slots
// and carry no group information, and some differ only in axis order.
[[nodiscard]] int site_symmetry_order(std::string_view symbol) {
  std::string key;
  std::ranges::copy(symbol |
                        std::views::filter([](char c) { return c != '.'; }),
                    std::back_inserter(key));
  static std::map<std::string_view, std::string_view> const kPermuted{
      {"2mm", "mm2"}, {"m2m", "mm2"}, {"-4m2", "-42m"}, {"-6m2", "-62m"}};
  if (auto const it = kPermuted.find(key); it != kPermuted.end()) {
    key = it->second;
  }
  static std::map<std::string_view, int> const kOrder{
      {"1", 1},     {"-1", 2},     {"2", 2},    {"m", 2},    {"2/m", 4},
      {"222", 4},   {"mm2", 4},    {"mmm", 8},  {"4", 4},    {"-4", 4},
      {"4/m", 8},   {"422", 8},    {"4mm", 8},  {"-42m", 8}, {"4/mmm", 16},
      {"3", 3},     {"-3", 6},     {"32", 6},   {"3m", 6},   {"-3m", 12},
      {"6", 6},     {"-6", 6},     {"6/m", 12}, {"622", 12}, {"6mm", 12},
      {"-62m", 12}, {"6/mmm", 24}, {"23", 12},  {"m-3", 24}, {"432", 24},
      {"-43m", 24}, {"m-3m", 48}};
  auto const it = kOrder.find(key);
  REQUIRE(it != kOrder.end());
  return it->second;
}
} // namespace

TEST_CASE("orbit-stabilizer invariant holds for all 230 space groups",
          "[group]") {
  // Three independently sourced quantities in one identity: tabulated
  // multiplicity, the stabiliser at a generic point of the tabulated locus, and
  // the orbit expanded from a generic seed. The site-symmetry symbol pins the
  // stabiliser's order a second way, so a pair mistabulated together is caught.
  Vector3d const seed{0.1234, 0.2718, 0.3142};
  int checked = 0;
  for (int number = 1; number <= 230; ++number) {
    auto const *sg =
        must(group::SpaceGroup::from_number(GroupFamily::space, number));
    auto const nops = static_cast<int>(sg->operations().size());
    REQUIRE(nops > 0);
    for (auto const &wp : sg->wyckoffs()) {
      INFO("space group " << number << ", position " << wp.letter() << " ("
                          << wp.site_symmetry() << ")");
      auto const order = static_cast<int>(wp.operations().size());
      // multiplicity * |site-symmetry group| == |conventional operations|
      REQUIRE(wp.multiplicity() * order == nops);
      // The expanded orbit of a generic point has exactly that many points.
      REQUIRE(wp.orbit(seed).rows() == wp.multiplicity());
      // The tabulated symbol names a point group of the stabiliser's order.
      REQUIRE(site_symmetry_order(wp.site_symmetry()) == order);
      ++checked;
    }
    // The general position is last, fully free, with trivial site symmetry.
    auto const &general = sg->wyckoffs().back();
    REQUIRE(general.degrees_of_freedom() == 3);
    REQUIRE(general.operations().size() == 1);
  }
  REQUIRE(checked > 1000);
}

TEST_CASE("SpaceGroup Pm-3m (221) matches ITA reference", "[group]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 221));
  REQUIRE(sg->number() == 221);
  REQUIRE(sg->operations().size() == 48);

  auto a = must(sg->wyckoff('a'));
  REQUIRE(a->multiplicity() == 1);
  REQUIRE(a->degrees_of_freedom() == 0);
  REQUIRE(a->site_symmetry() == "m-3m");

  // 1a has no free coordinate: its orbit is the single fixed point, and every
  // point projects onto it whatever the free parameters say.
  Positions const orbit = a->orbit(Vector3d{0.3, 0.4, 0.5});
  REQUIRE(orbit.rows() == 1);
  REQUIRE(a->canonical(Vector3d{0.3, 0.4, 0.5}).isApprox(a->sample({})));

  REQUIRE(sg->wyckoffs().back().multiplicity() == 48);
}

TEST_CASE("Wyckoff orbit expansion respects multiplicity", "[group]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 225)); // Fm-3m
  auto const &general = sg->wyckoffs().back();
  Positions const orbit = general.orbit(Vector3d{0.11, 0.23, 0.37});
  REQUIRE(orbit.rows() == general.multiplicity());
}

TEST_CASE("from_number rejects out-of-range numbers", "[group]") {
  REQUIRE(errored(
      [] { return group::SpaceGroup::from_number(GroupFamily::space, 231); }));
}

TEST_CASE("crystal generation round-trips through the analyzer", "[generate]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 225)); // Fm-3m
  // NaCl: Na (type 11) and Cl (type 17), four of each (the 4a / 4b orbits).
  generate::Composition const comp{{11, 4}, {17, 4}};
  REQUIRE(generate::Generator{*sg}.compatible(comp));

  auto gen = must(generate::Generator{*sg, {.seed = 42u}}(comp));
  REQUIRE(gen.cell.size() == 8);
  std::set<int> kinds(gen.cell.types().begin(), gen.cell.types().end());
  REQUIRE(kinds == std::set<int>{11, 17});

  // The assembled structure actually has the requested symmetry.
  auto analyzer = analysis::SymmetryAnalyzer::from_cell(gen.cell);
  REQUIRE(data::spacegroup_type(must(analyzer.hall())).number == 225);
}

TEST_CASE("generated structures recover their target space group",
          "[generate]") {
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
    auto const *sg =
        must(group::SpaceGroup::from_number(GroupFamily::space, c.number));
    auto gen = must(generate::Generator{*sg, {.seed = 7u}}(c.comp));
    auto analyzer = analysis::SymmetryAnalyzer::from_cell(gen.cell);
    REQUIRE(data::spacegroup_type(must(analyzer.hall())).number == c.number);
  }
}

TEST_CASE("generated structures are free of interatomic clashes",
          "[generate]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 225));
  generate::Composition const comp{{11, 4}, {17, 4}};
  auto gen = must(generate::Generator{*sg, {.seed = 99u}}(comp));
  // The builder only returns distance-valid cells; confirm independently.
  REQUIRE(generate::distances_valid(gen.cell));
}

TEST_CASE("generation is deterministic in the seed", "[generate]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 225));
  generate::Composition const comp{{11, 4}, {17, 4}};
  generate::Generator const gen{*sg, {.seed = 123u}};
  auto a = must(gen(comp));
  auto b = must(gen(comp));
  REQUIRE(a.cell.positions().isApprox(b.cell.positions()));
  REQUIRE(a.cell.lattice().matrix().isApprox(b.cell.lattice().matrix()));
}

TEST_CASE("incompatible composition is rejected", "[generate]") {
  auto const *sg = must(group::SpaceGroup::from_number(
      GroupFamily::space, 225)); // Fm-3m, smallest mult 4
  // Three atoms cannot fill any combination of multiplicity-4-or-more orbits.
  REQUIRE_FALSE(
      generate::Generator{*sg}.compatible(generate::Composition{{1, 3}}));
}

TEST_CASE("fixed sites pin their Wyckoff letters and are counted",
          "[generate]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 225));
  generate::Composition const comp{{11, 4}, {17, 4}};
  generate::GenerateOptions options{.seed = 5u};
  options.sites = {{.type = 11, .letter = 'a'}, {.type = 17, .letter = 'b'}};
  generate::Generator const gen{*sg, options};

  // Everything is pinned: exactly one assignment, the fixed one.
  auto const all = gen.assignments(comp);
  REQUIRE(all.size() == 1);
  CHECK(all.front()[0].position->letter() == 'a');
  CHECK(all.front()[1].position->letter() == 'b');

  auto const generated = must(gen(comp));
  REQUIRE(generated.assignment.size() == 2);
  CHECK(generated.assignment[0].coordinate.has_value());
  CHECK(data::spacegroup_type(must(test::dataset_of(generated.cell)).hall)
            .number == 225);

  // A letter the group lacks is an error; overshooting the composition is
  // an incompatibility.
  options.sites = {{.type = 11, .letter = 'z'}};
  CHECK(test::errored([&] { return generate::Generator{*sg, options}(comp); }));
  options.sites = {{.type = 11, .letter = 'a'}, {.type = 11, .letter = 'b'}};
  CHECK_FALSE(generate::Generator{*sg, options}.compatible(comp));
}

TEST_CASE("a fixed coordinate is kept, projected onto its locus",
          "[generate]") {
  auto const *sg = must(group::SpaceGroup::from_number(GroupFamily::space, 1));
  Vector3d const pinned(0.1, 0.2, 0.3);
  generate::GenerateOptions options{.seed = 3u};
  options.sites = {{.type = 6, .letter = 'a', .coordinate = pinned}};
  auto const generated =
      must(generate::Generator{*sg, options}(generate::Composition{{6, 2}}));
  REQUIRE(generated.cell.size() == 2);
  CHECK(generated.assignment.front().coordinate->isApprox(pinned));
  CHECK(generated.cell.position(0).isApprox(pinned));
}

TEST_CASE("a caller's lattice is used as given or rejected", "[generate]") {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, 221));
  generate::Composition const comp{{55, 1}, {17, 1}}; // CsCl
  generate::GenerateOptions options{.seed = 1u};
  options.lattice = Lattice{Matrix3d::Identity() * 4.1};
  auto const generated = must(generate::Generator{*sg, options}(comp));
  CHECK(generated.cell.lattice().matrix().isApprox(Matrix3d::Identity() * 4.1));

  Matrix3d sheared = Matrix3d::Identity() * 4.1;
  sheared(0, 1) = 0.7; // no longer cubic
  options.lattice = Lattice{sheared};
  CHECK(test::errored([&] { return generate::Generator{*sg, options}(comp); }));
}

TEST_CASE("DistanceTolerance presets and overrides", "[generate]") {
  using generate::DistanceTolerance;
  DistanceTolerance const covalent;
  auto const metallic = DistanceTolerance::preset(data::RadiusKind::metallic);
  CHECK(covalent.min_distance(26, 26) ==
        Catch::Approx(2 * 0.7 * *data::covalent_radius(26)));
  CHECK(metallic.min_distance(26, 26) ==
        Catch::Approx(2 * 0.7 * *data::radius(data::RadiusKind::metallic, 26)));
  CHECK(metallic.min_distance(26, 26) != covalent.min_distance(26, 26));
  // He has no metallic radius: the fallback stands in.
  CHECK(metallic.radius(2) == 1.0);
  CHECK(covalent.radius(2) == *data::covalent_radius(2));

  DistanceTolerance pinned;
  pinned.set(17, 11, 3.0);
  CHECK(pinned.min_distance(11, 17) == 3.0);
  CHECK(pinned.min_distance(17, 11) == 3.0);
  CHECK(pinned.min_distance(11, 11) == covalent.min_distance(11, 11));
}

TEST_CASE("orbit-stabilizer invariant holds for all 80 layer groups",
          "[layergen]") {
  // The same three-way check as for the space groups; the orbit folds only
  // the periodic plane, so a c-flipped image stays at -z rather than 1-z.
  Vector3d const seed{0.1234, 0.2718, 0.3142};
  for (int number = 1; number <= 80; ++number) {
    auto const *lg =
        must(group::SpaceGroup::from_number(GroupFamily::layer, number));
    REQUIRE(lg->number() == number);
    REQUIRE(lg->hall().family() == GroupFamily::layer);
    auto const nops = static_cast<int>(lg->operations().size());
    REQUIRE(nops > 0);
    for (auto const &wp : lg->wyckoffs()) {
      INFO("layer group " << number << ", position " << wp.letter() << " ("
                          << wp.site_symmetry() << ")");
      auto const order = static_cast<int>(wp.operations().size());
      REQUIRE(wp.multiplicity() * order == nops);
      REQUIRE(wp.orbit(seed, aperiodic_along(2)).rows() == wp.multiplicity());
      REQUIRE(site_symmetry_order(wp.site_symmetry()) == order);
    }
  }
}

TEST_CASE("generated layer structures carry their full layer symmetry",
          "[layergen]") {
  // A generated 2D-periodic structure must be invariant under EVERY operation
  // of its layer group, under aperiodic-aware overlap. Covers the c-flipping
  // groups, whose orbits must build without folding the non-periodic axis.
  for (int number : {1, 2, 19, 49, 61, 65}) {
    INFO("layer group " << number);
    auto const *lg =
        must(group::SpaceGroup::from_number(GroupFamily::layer, number));
    int const m = lg->wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    // p4/mmm (61) puts 32 atoms in a thin slab, close enough to its packing
    // limit that a 50-attempt budget clears it only about three times in four;
    // the budget, not the metric, is what this test needs raised.
    auto gen = must(generate::Generator{
        *lg,
        {.scale = 4.0,
         .seed = 13u,
         .attempts_per_combination = 400,
         .placement = generate::Placement::general_only}}(comp));
    REQUIRE(aperiodic_axis(gen.cell.periodicity()) == 2);
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(generate::distances_valid(gen.cell));

    OverlapChecker checker(gen.cell, 1e-3);
    for (auto const &op : lg->operations()) {
      REQUIRE(checker.check_total_overlap(op.translation, op.rotation));
    }
  }
}

TEST_CASE("layer crystal round-trips through the layer dataset", "[layergen]") {
  // End-to-end self-validation through the layer determination. Restricted to
  // c-preserving layer groups (one atomic plane): the layer path currently
  // does not resolve multi-level (c-flipping) layers, so those are covered by
  // the direct op-invariance test above rather than a determination round-trip.
  for (int number : {1, 49, 55, 65}) {
    INFO("layer group " << number);
    auto const *lg =
        must(group::SpaceGroup::from_number(GroupFamily::layer, number));
    int const m = lg->wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::Generator{
        *lg,
        {.scale = 4.0,
         .seed = 13u,
         .attempts_per_combination = 400,
         .placement = generate::Placement::general_only}}(comp));
    auto ds = must(test::dataset_of(
        gen.cell.with_periodicity(aperiodic_along(2)), {1e-4}));
    REQUIRE(data::spacegroup_type(ds.hall).number ==
            number); // exact layer-group recovery
  }
}

TEST_CASE("SymmetryAnalyzer memoizes a consistent dataset", "[analysis]") {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  Cell const cell{Lattice{Matrix3d::Identity()}, pos, Types{1}};

  auto analyzer = analysis::SymmetryAnalyzer::from_cell(cell);
  int const first = data::spacegroup_type(must(analyzer.hall())).number;
  int const cached = data::spacegroup_type(must(analyzer.hall())).number;
  REQUIRE(first == 221); // simple cubic, one atom -> Pm-3m
  REQUIRE(first == cached);

  // Projections agree with the dataset.
  auto ds = must(analyzer.dataset());
  REQUIRE(data::spacegroup_type(ds.hall).number == first);
  REQUIRE(must(analyzer.operations()).size() == ds.operations.size());
}

namespace {
// A structure carries its group's full symmetry iff every operation maps the
// atom set onto itself type-for-type, under the cell's OWN periodicity: a flip
// along an aperiodic axis lands at -a, never the wrapped 1 - a. That covers the
// cluster (nothing folded) and the rod (only its periodic axis).
bool is_invariant(Cell const &cell, std::span<SymmetryOperation const> ops,
                  double tol = 1e-4) {
  Index const n = cell.size();
  for (auto const &op : ops) {
    for (Index i = 0; i < n; ++i) {
      Vector3d const image = op.apply(cell.position(i));
      bool matched = false;
      for (Index j = 0; j < n && !matched; ++j) {
        matched = cell.type(i) == cell.type(j) &&
                  minimal_image(Vector3d(image - cell.position(j)),
                                cell.periodicity())
                          .cwiseAbs()
                          .maxCoeff() < tol;
      }
      if (!matched) {
        return false;
      }
    }
  }
  return true;
}
} // namespace

TEST_CASE("orbit-stabilizer invariant holds for all 32 point groups",
          "[cluster]") {
  for (int number = 1; number <= 32; ++number) {
    INFO("point group " << number);
    auto pg = must(group::PointGroup::from_number(number));
    REQUIRE(pg.number() == number);
    int const order = pg.order();
    REQUIRE(order > 0);
    REQUIRE_FALSE(pg.wyckoffs().empty());

    // The origin is always a fixed point (multiplicity 1); the general position
    // (last) has multiplicity equal to the group order.
    REQUIRE(pg.wyckoffs().front().multiplicity() == 1);
    REQUIRE(pg.wyckoffs().back().multiplicity() == order);
    REQUIRE(pg.wyckoffs().back().degrees_of_freedom() == 3);

    for (auto const &wp : pg.wyckoffs()) {
      REQUIRE(wp.multiplicity() * static_cast<int>(wp.operations().size()) ==
              order);
    }
  }
}

TEST_CASE("point-group orders match the textbook values", "[cluster]") {
  auto order_of = [](int number) {
    return must(group::PointGroup::from_number(number)).order();
  };
  REQUIRE(order_of(1) == 1);   // 1   (C1)
  REQUIRE(order_of(2) == 2);   // -1  (Ci)
  REQUIRE(order_of(5) == 4);   // 2/m (C2h)
  REQUIRE(order_of(8) == 8);   // mmm (D2h)
  REQUIRE(order_of(16) == 3);  // 3   (C3)
  REQUIRE(order_of(20) == 12); // -3m (D3d)
  REQUIRE(order_of(25) == 12); // 6mm (C6v)
  REQUIRE(order_of(27) == 24); // 6/mmm (D6h)
  REQUIRE(order_of(32) == 48); // m-3m (Oh)
  auto const symbol_of = [](int number) {
    return must(group::PointGroup::from_number(number)).schoenflies();
  };
  REQUIRE(symbol_of(5) == "C2h");
  REQUIRE(symbol_of(32) == "Oh");
}

TEST_CASE("from_number rejects out-of-range point groups", "[cluster]") {
  REQUIRE(errored([] { return group::PointGroup::from_number(33); }));
}

TEST_CASE("generated clusters carry their full point-group symmetry",
          "[cluster]") {
  // The direct self-validation: a generated 0D cluster must be invariant under
  // EVERY operation of its point group, in Cartesian space. Spans all crystal
  // systems including the trigonal/hexagonal metrics, where the integer
  // operations are only isometries in the correct (non-orthogonal) basis.
  for (int number : {1, 2, 8, 16, 20, 25, 32}) {
    INFO("point group " << number);
    auto pg = must(group::PointGroup::from_number(number));
    int const m = pg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::Generator{
        pg,
        {.scale = 3.0,
         .seed = 17u,
         .placement = generate::Placement::general_only}}(comp));
    REQUIRE(gen.cell.periodicity() == none_periodic());
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(generate::distances_valid(gen.cell));
    REQUIRE(is_invariant(gen.cell, pg.operations(), 1e-6));
  }
}

TEST_CASE("cluster generation is deterministic in the seed", "[cluster]") {
  auto pg = must(group::PointGroup::from_number(32)); // m-3m
  generate::Composition const comp{{6, 48}};
  generate::Generator const gen{pg, {.seed = 123u}};
  auto a = must(gen(comp));
  auto b = must(gen(comp));
  REQUIRE(a.cell.positions().isApprox(b.cell.positions()));
  REQUIRE(a.cell.types() == b.cell.types());
}

TEST_CASE("incompatible cluster composition is rejected", "[cluster]") {
  auto pg = must(group::PointGroup::from_number(32)); // m-3m, order 48
  // 5 atoms cannot tile the available multiplicities (1, ..., 48) of m-3m.
  generate::Composition const comp{{6, 5}};
  REQUIRE(errored([&] { return generate::Generator{pg}(comp); }));
}

TEST_CASE("orbit-stabilizer invariant holds for all 75 rod groups", "[rod]") {
  REQUIRE(data::num_rod_groups() == 75);
  for (int number = 1; number <= 75; ++number) {
    INFO("rod group " << number);
    auto rg = must(group::RodGroup::from_number(number));
    REQUIRE(rg.number() == number);
    REQUIRE(rg.periodic_axis() == 2);
    int const order = rg.order();
    REQUIRE(order > 0);
    REQUIRE_FALSE(rg.wyckoffs().empty());

    // The general position (last) is fully free with multiplicity == order;
    // every position satisfies multiplicity * |site symmetry| == order.
    REQUIRE(rg.wyckoffs().back().multiplicity() == order);
    REQUIRE(rg.wyckoffs().back().degrees_of_freedom() == 3);
    for (auto const &wp : rg.wyckoffs()) {
      REQUIRE(wp.multiplicity() * static_cast<int>(wp.operations().size()) ==
              order);
    }
  }
}

TEST_CASE("generated rod structures carry their full rod symmetry", "[rod]") {
  // Direct self-validation of the general-position path: the generated 1D
  // structure must be invariant under EVERY rod operation, folding only the
  // periodic axis. Spans the trivial group, a perpendicular 2-fold
  // (c-flipping), a screw axis, and the higher-symmetry axial groups.
  for (int number : {1, 3, 13, 23, 24, 53}) {
    INFO("rod group " << number);
    auto rg = must(group::RodGroup::from_number(number));
    int const m = rg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::Generator{
        rg,
        {.scale = 2.0,
         .seed = 31u,
         .placement = generate::Placement::general_only}}(comp));
    REQUIRE(gen.cell.periodicity() == periodic_along(2));
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(is_invariant(gen.cell, rg.operations()));
  }
}

TEST_CASE("rod special positions are derived and generate invariant structures",
          "[rod]") {
  // The affine fixed-locus arrangement yields special positions (atoms on the
  // rod axis / cross-section special points), not just the general one. Placing
  // atoms on the most special position must still produce a structure invariant
  // under every operation.
  for (int number : {8, 13, 23, 37}) { // p112, p222, p4, p-42m
    auto rg = must(group::RodGroup::from_number(number));
    INFO("rod group " << number << " (" << rg.symbol() << ") with "
                      << rg.wyckoffs().size() << " Wyckoff positions");
    REQUIRE(rg.wyckoffs().size() >= 2); // at least one special + the general
    auto const &special = rg.wyckoffs().front();
    REQUIRE(special.multiplicity() < rg.order()); // genuinely special

    generate::Composition const comp{{6, special.multiplicity()}};
    auto gen = must(generate::Generator{rg, {.scale = 2.0, .seed = 5u}}(comp));
    REQUIRE(is_invariant(gen.cell, rg.operations()));
  }
}

TEST_CASE("rod generation is deterministic in the seed", "[rod]") {
  auto rg = must(group::RodGroup::from_number(23)); // p4
  int const m = rg.wyckoffs().back().multiplicity();
  generate::Composition const comp{{6, m}};
  generate::Generator const gen{rg, {.seed = 7u}};
  auto a = must(gen(comp));
  auto b = must(gen(comp));
  REQUIRE(a.cell.positions().isApprox(b.cell.positions()));
  REQUIRE(a.cell.types() == b.cell.types());
}

TEST_CASE("warmup builds the shared per-setting caches", "[warmup]") {
  // Priming is idempotent and races the ordinary first-use path safely: every
  // caller ends up with the same flyweight whichever built it.
  warmup_async(Warm::layer_groups).get();
  warmup(Warm::all);
  for (auto const family : {GroupFamily::space, GroupFamily::layer}) {
    for (int index : {1, hall_settings(family)}) {
      auto const key = *HallNumber::of(family, index);
      auto const &group = group::SpaceGroup::of(key);
      CHECK(&group == &group::SpaceGroup::of(key)); // shared, not rebuilt
      CHECK_FALSE(group.wyckoffs().empty());
    }
  }
}
