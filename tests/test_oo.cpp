// Tests for the object-oriented layer: the analysis facades (SymmetryAnalyzer),
// the standalone group objects (SpaceGroup / WyckoffPosition), and crystal
// generation. The group tests assert the orbit-stabilizer invariant across all
// 230 space groups; the generation test round-trips a generated cell back
// through the analyzer.
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include "core/overlap.hpp"
#include "data/rod_database.hpp"
#include <cppcrystal/dataset.hpp>
#include <cppcrystal/generate/crystal_builder.hpp>
#include <cppcrystal/generate/distance_check.hpp>
#include <cppcrystal/generate/rod_crystal.hpp>
#include <cppcrystal/group/point_group.hpp>
#include <cppcrystal/group/rod_group.hpp>
#include <cppcrystal/group/space_group.hpp>
#include <cppcrystal/group/subgroup_graph.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

using namespace cppcrystal;

namespace {
// Unwrap a move-only Result<T>, failing the test on an error result.
using cppcrystal::test::errored;
using cppcrystal::test::must;
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
  REQUIRE(errored([] { return group::SpaceGroup::from_number(231); }));
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
  REQUIRE(a.cell.lattice().matrix().isApprox(b.cell.lattice().matrix()));
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
        {.scale = 4.0, .seed = 13u, .general_position_only = true}));
    REQUIRE(aperiodic_axis(gen.cell.periodicity()) == 2);
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(generate::distances_valid(gen.cell));

    OverlapChecker checker(gen.cell, 1e-3);
    for (auto const &op : lg.operations()) {
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
    auto lg = must(group::SpaceGroup::from_layer_number(number));
    int const m = lg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::random_layer_crystal(
        lg, comp,
        {.scale = 4.0, .seed = 13u, .general_position_only = true}));
    auto ds = must(get_dataset(gen.cell.with_periodicity(aperiodic_along(2)),
                               {1e-4}));
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
  Cell const cell{Lattice{Matrix3d::Identity()}, pos, Types{1}};

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

namespace {
// A generated cluster carries the full point-group symmetry iff, in Cartesian
// space, every operation R_cart = basis . R . basis^-1 maps the atom set onto
// itself (each atom onto an atom of the same type).
bool cluster_is_invariant(generate::GeneratedCluster const &gen,
                          group::PointGroup const &pg, double tol = 1e-6) {
  Matrix3d const inv = gen.basis.inverse();
  Index const n = static_cast<Index>(gen.types.size());
  for (auto const &op : pg.operations()) {
    Matrix3d const rc = gen.basis * op.rotation.cast<double>() * inv;
    for (Index i = 0; i < n; ++i) {
      Vector3d const image = rc * gen.coordinates.row(i).transpose();
      bool matched = false;
      for (Index j = 0; j < n && !matched; ++j) {
        if (gen.types[static_cast<std::size_t>(i)] ==
                gen.types[static_cast<std::size_t>(j)] &&
            (image - gen.coordinates.row(j).transpose()).norm() < tol) {
          matched = true;
        }
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
      REQUIRE(wp.multiplicity() *
                  static_cast<int>(wp.operations().size()) ==
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

    auto gen = must(generate::random_cluster(
        pg, comp,
        {.scale = 3.0, .seed = 17u, .general_position_only = true}));
    REQUIRE(gen.types.size() == static_cast<std::size_t>(2 * m));
    REQUIRE(generate::cluster_distances_valid(gen.coordinates, gen.types));
    REQUIRE(cluster_is_invariant(gen, pg));
  }
}

TEST_CASE("cluster generation is deterministic in the seed", "[cluster]") {
  auto pg = must(group::PointGroup::from_number(32)); // m-3m
  generate::Composition const comp{{6, 48}};
  auto a = must(generate::random_cluster(pg, comp, {.seed = 123u}));
  auto b = must(generate::random_cluster(pg, comp, {.seed = 123u}));
  REQUIRE(a.coordinates.isApprox(b.coordinates));
  REQUIRE(a.types == b.types);
}

TEST_CASE("incompatible cluster composition is rejected", "[cluster]") {
  auto pg = must(group::PointGroup::from_number(32)); // m-3m, order 48
  // 5 atoms cannot tile the available multiplicities (1, ..., 48) of m-3m.
  generate::Composition const comp{{6, 5}};
  REQUIRE(errored([&] { return generate::random_cluster(pg, comp); }));
}

namespace {
// A generated rod carries the full rod-group symmetry iff every operation maps
// the fractional atom set onto itself, folding ONLY the periodic axis (the two
// aperiodic axes are compared raw — a flip a -> -a must land at the Cartesian
// image, not the wrapped 1 - a).
bool rod_is_invariant(generate::GeneratedRodCrystal const &gen,
                      group::RodGroup const &rg, double tol = 1e-4) {
  int const axis = rg.periodic_axis();
  Cell const &cell = gen.cell;
  Index const n = cell.size();
  for (auto const &op : rg.operations()) {
    for (Index i = 0; i < n; ++i) {
      Vector3d const image = op.apply(cell.position(i));
      bool matched = false;
      for (Index j = 0; j < n && !matched; ++j) {
        if (cell.type(i) != cell.type(j)) {
          continue;
        }
        Vector3d d = image - cell.position(j);
        d[axis] -= std::round(d[axis]); // fold only the periodic axis
        if (d.cwiseAbs().maxCoeff() < tol) {
          matched = true;
        }
      }
      if (!matched) {
        return false;
      }
    }
  }
  return true;
}
} // namespace

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
      REQUIRE(wp.multiplicity() *
                  static_cast<int>(wp.operations().size()) ==
              order);
    }
  }
}

TEST_CASE("generated rod structures carry their full rod symmetry", "[rod]") {
  // Direct self-validation of the general-position path: the generated 1D
  // structure must be invariant under EVERY rod operation, folding only the
  // periodic axis. Spans the trivial group, a perpendicular 2-fold (c-flipping),
  // a screw axis, and the higher-symmetry axial groups.
  for (int number : {1, 3, 13, 23, 24, 53}) {
    INFO("rod group " << number);
    auto rg = must(group::RodGroup::from_number(number));
    int const m = rg.wyckoffs().back().multiplicity();
    generate::Composition const comp{{6, m}, {7, m}};

    auto gen = must(generate::random_rod_crystal(
        rg, comp,
        {.scale = 2.0, .seed = 31u, .general_position_only = true}));
    REQUIRE(gen.cell.periodicity() == CellPeriodicity{AxisKind::aperiodic,
                                                      AxisKind::aperiodic,
                                                      AxisKind::periodic});
    REQUIRE(gen.cell.size() == static_cast<Index>(2 * m));
    REQUIRE(rod_is_invariant(gen, rg));
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
    auto gen = must(generate::random_rod_crystal(
        rg, comp, {.scale = 2.0, .seed = 5u}));
    REQUIRE(rod_is_invariant(gen, rg));
  }
}

TEST_CASE("rod generation is deterministic in the seed", "[rod]") {
  auto rg = must(group::RodGroup::from_number(23)); // p4
  int const m = rg.wyckoffs().back().multiplicity();
  generate::Composition const comp{{6, m}};
  auto a = must(generate::random_rod_crystal(rg, comp, {.seed = 7u}));
  auto b = must(generate::random_rod_crystal(rg, comp, {.seed = 7u}));
  REQUIRE(a.cell.positions().isApprox(b.cell.positions()));
  REQUIRE(a.cell.types() == b.cell.types());
}
