// The maximal-subgroup graph: table invariants, textbook relations, the BGL
// face. Checks needing the derived (klassengleiche) table SKIP without it.

#include <seitz/data/spg_database.hpp>
#include <seitz/group/space_group.hpp>
#include <seitz/group/subgroup_graph.hpp>

#include "core/matrix_order.hpp"
#include "math/integer_matrix.hpp"

#include "helpers.hpp"

#include <boost/graph/graph_traits.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <iterator>
#include <ranges>
#include <vector>

using namespace seitz;
using group::SubgroupEdge;
using group::SubgroupGraph;
using group::SubgroupKind;
using test::must;

namespace {

// Point-group order of a space group: distinct rotation parts of its
// operations.
int pg_order(int number) {
  auto const *sg =
      must(group::SpaceGroup::from_number(GroupFamily::space, number));
  return static_cast<int>(
      rotation_set(sg->operations(), &SymmetryOperation::rotation).size());
}

// Lattice points per conventional cell: the pure translations among the
// setting's operations.
int lattice_points(HallNumber hall) {
  return static_cast<int>(std::ranges::count_if(
      group::SpaceGroup::of(hall).operations(),
      [](SymmetryOperation const &op) { return op.rotation.isIdentity(); }));
}

std::vector<SubgroupEdge> edges_of(std::ranges::input_range auto ids) {
  return {std::from_range, ids | std::views::transform(SubgroupGraph::edge)};
}

bool has_relation(std::vector<SubgroupEdge> const &edges, int sub, int index) {
  return std::ranges::any_of(edges, [&](SubgroupEdge const &e) {
    return e.sub == sub && e.index == index;
  });
}

// Whether the table is the derived one (klassengleiche edges present) rather
// than the translationengleiche-only bootstrap with placeholder bases.
bool derived_table() {
  auto const ids =
      std::views::iota(0, static_cast<int>(SubgroupGraph::num_edges()));
  return std::ranges::any_of(ids, [](int id) {
    return SubgroupGraph::edge(id).kind == SubgroupKind::klassengleiche;
  });
}

constexpr char const *kNotDerived =
    "the klassengleiche table has not been derived yet: run "
    "tools/generate_subgroup_relations.py and rebuild";

// Whether every operation of the subgroup's setting, taken to the parent frame
// (x_G = P x_H + s), is an operation of the parent: strictly H <= G.
bool subgroup_is_in_parent(SubgroupEdge const &e) {
  auto const &parent =
      group::SpaceGroup::of(*data::default_hall<GroupFamily::space>(e.super))
          .operations();
  Matrix3d const p_inv = e.basis.inverse();
  auto const in_parent = [&](SymmetryOperation const &op) {
    Matrix3d const rotation = e.basis * op.rotation.cast<double>() * p_inv;
    Vector3d const translation =
        e.basis * op.translation + e.origin - rotation * e.origin;
    return math::is_int_matrix(rotation, 1e-6) &&
           std::ranges::any_of(parent, [&](SymmetryOperation const &d) {
             return d.rotation == math::round_to_int(rotation) &&
                    math::nearest_offset(
                        Vector3d(translation - d.translation))
                            .norm() < 1e-6;
           });
  };
  return std::ranges::all_of(group::SpaceGroup::of(e.hall).operations(),
                             in_parent);
}

} // namespace

TEST_CASE("Pm-3m maximal t-subgroups match the textbook relations",
          "[subgroup]") {
  auto const subs = edges_of(SubgroupGraph::maximal_subgroups(
      221, SubgroupKind::translationengleiche)); // Oh
  // The five maximal translationengleiche subgroups of Pm-3m (ITA): Pm-3 (200,
  // Th), P432 (207, O), P-43m (215, Td) at index 2; P4/mmm (123, D4h) at index
  // 3; R-3m (166, D3d, along a body diagonal) at index 4 -- maximal subgroups
  // of a non-nilpotent group may have composite index.
  REQUIRE(has_relation(subs, 200, 2));
  REQUIRE(has_relation(subs, 207, 2));
  REQUIRE(has_relation(subs, 215, 2));
  REQUIRE(has_relation(subs, 123, 3));
  REQUIRE(has_relation(subs, 166, 4));
  REQUIRE(subs.size() == 5);

  if (derived_table()) {
    // R-3m's hexagonal cell holds three cubic lattice points.
    auto const r3m = std::ranges::find(subs, 166, &SubgroupEdge::sub);
    REQUIRE(math::nint(r3m->basis.determinant() *
                       lattice_points(
                           *data::default_hall<GroupFamily::space>(221))) ==
            lattice_points(r3m->hall));
  }
}

TEST_CASE("P-1 has the single maximal t-subgroup P1", "[subgroup]") {
  auto const subs = edges_of(
      SubgroupGraph::maximal_subgroups(2, SubgroupKind::translationengleiche));
  REQUIRE(subs.size() == 1);
  REQUIRE(subs.front().sub == 1);
  REQUIRE(subs.front().index == 2);
  // P1 has no translationengleiche subgroup at all.
  REQUIRE(std::ranges::empty(SubgroupGraph::maximal_subgroups(
      1, SubgroupKind::translationengleiche)));
}

TEST_CASE("every edge is consistent with its endpoints", "[subgroup]") {
  bool const derived = derived_table();
  int edges = 0;
  for (int n = 1; n <= group::kNumSpaceGroups; ++n) {
    for (SubgroupEdge const &e : edges_of(SubgroupGraph::maximal_subgroups(n))) {
      INFO("edge " << e.id << ": " << n << " -> " << e.sub << " (hall "
                   << e.hall.index() << ")");
      REQUIRE(e.super == n);
      REQUIRE(e.index >= 2);
      REQUIRE(e.basis.determinant() > 0);
      if (e.kind == SubgroupKind::translationengleiche) {
        // The lattice is kept, so the index is the point-group order ratio.
        REQUIRE(e.sub != n);
        REQUIRE(pg_order(n) == e.index * pg_order(e.sub));
      } else {
        // The point group is kept; only a small sublattice is dropped.
        REQUIRE(pg_order(n) == pg_order(e.sub));
        REQUIRE((e.index >= 2 && e.index <= 4));
      }
      if (derived) {
        // Cell volumes: |det P| . m_G == (index) . m_H, the index being 1 for
        // a translationengleiche edge (same lattice).
        int const lattice_index =
            e.kind == SubgroupKind::klassengleiche ? e.index : 1;
        REQUIRE(math::nint(e.basis.determinant() *
                           lattice_points(
                               *data::default_hall<GroupFamily::space>(n))) ==
                lattice_index * lattice_points(e.hall));
        // The setting is a genuine embedding: every operation of the subgroup
        // is an operation of the parent.
        REQUIRE(subgroup_is_in_parent(e));
      }
      // The in-edges are the same edge set read the other way round.
      REQUIRE(std::ranges::contains(SubgroupGraph::minimal_supergroups(e.sub),
                                    e.id));
      ++edges;
    }
  }
  REQUIRE(edges > 200); // the graph is richly connected
}

TEST_CASE("known klassengleiche relations are present", "[subgroup]") {
  if (!derived_table()) {
    SKIP(kNotDerived);
  }
  auto const k = [](int n) {
    return edges_of(
        SubgroupGraph::maximal_subgroups(n, SubgroupKind::klassengleiche));
  };
  // Pm-3m -> Fm-3m with the cell doubled along every axis.
  auto const from_221 = k(221);
  auto const fm3m = std::ranges::find(from_221, 225, &SubgroupEdge::sub);
  REQUIRE(fm3m != from_221.end());
  REQUIRE(fm3m->index == 2);
  // Doubled along every axis, up to the cubic setting the determination chose:
  // |basis| / 2 is a permutation matrix.
  Matrix3d const halved = fm3m->basis.cwiseAbs() / 2;
  REQUIRE(halved.colwise().sum().isApprox(Eigen::RowVector3d::Ones()));
  REQUIRE(halved.rowwise().sum().isApprox(Vector3d::Ones()));
  REQUIRE(halved.cwiseProduct(halved).isApprox(halved));
  // P-31m -> R-3m and P-3 -> R-3 (the R lattice as an index-3 sublattice of
  // the hexagonal one), Pmmm -> Cmmm at index 2.
  REQUIRE(has_relation(k(162), 166, 3));
  REQUIRE(has_relation(k(147), 148, 3));
  REQUIRE(has_relation(k(47), 65, 2));
  // Isomorphic subgroups: P1 and P-1 in cells of index 2 and 3; an index-4
  // sublattice always sits inside an index-2 one, so no index-4 self-edge.
  for (int n : {1, 2}) {
    REQUIRE(has_relation(k(n), n, 2));
    REQUIRE(has_relation(k(n), n, 3));
    REQUIRE_FALSE(has_relation(k(n), n, 4));
  }
}

TEST_CASE("reachability and symmetry-breaking paths", "[subgroup]") {
  REQUIRE(SubgroupGraph::is_subgroup(221, 221)); // reflexive
  REQUIRE(SubgroupGraph::is_subgroup(1, 221));   // P1 lies under Pm-3m
  REQUIRE_FALSE(SubgroupGraph::is_subgroup(221, 1));
  REQUIRE_FALSE(SubgroupGraph::is_subgroup(0, 1));

  auto const chain = SubgroupGraph::path(221, 1);
  REQUIRE(chain.has_value());
  REQUIRE_FALSE(chain->empty());
  auto const edges = edges_of(*chain);
  REQUIRE(edges.front().super == 221);
  REQUIRE(edges.back().sub == 1);
  // Consecutive edges chain: each starts where the previous one ends.
  for (auto const &[a, b] : edges | std::views::pairwise) {
    REQUIRE(a.sub == b.super);
  }

  auto const t_only =
      SubgroupGraph::path(221, 1, SubgroupKind::translationengleiche);
  REQUIRE(t_only.has_value());
  REQUIRE(std::ranges::all_of(edges_of(*t_only), [](SubgroupEdge const &e) {
    return e.kind == SubgroupKind::translationengleiche;
  }));

  REQUIRE(SubgroupGraph::path(5, 5)->empty());
  REQUIRE_FALSE(SubgroupGraph::path(1, 221).has_value());
  REQUIRE_FALSE(SubgroupGraph::path(0, 1).has_value());

  if (derived_table()) {
    // Fm-3m keeps Pm-3m's point group, so it is one klassengleiche step away
    // and no translationengleiche chain reaches it.
    REQUIRE(SubgroupGraph::path(221, 225, SubgroupKind::klassengleiche)->size() ==
            1);
    REQUIRE_FALSE(
        SubgroupGraph::path(221, 225, SubgroupKind::translationengleiche)
            .has_value());
  }
}

TEST_CASE("the table is a Boost.Graph", "[subgroup]") {
  auto const &g = group::kSubgroupGraph;
  REQUIRE(num_vertices(g) == group::kNumSpaceGroups + 1);
  REQUIRE(num_edges(g) == SubgroupGraph::num_edges());

  auto const [first, last] = out_edges(221, g);
  REQUIRE(std::distance(first, last) ==
          std::ranges::distance(SubgroupGraph::maximal_subgroups(221)));
  for (int const e : std::ranges::subrange(first, last)) {
    REQUIRE(source(e, g) == 221);
    REQUIRE(target(e, g) == SubgroupGraph::edge(e).sub);
    REQUIRE(g[e].super == 221);
  }
  // The edge bundle's members are property maps.
  auto const index = get(&data::SubgroupRelation::index, g);
  REQUIRE(get(index, *first) == SubgroupGraph::edge(*first).index);
  // A kind view hides the other kind.
  group::SubgroupView const t_view{
      g, group::detail::KindIs{SubgroupKind::translationengleiche}};
  auto const [tf, tl] = out_edges(221, t_view);
  REQUIRE(std::distance(tf, tl) == 5);
  // In-edges: P1 is a maximal subgroup of every index-2 group at least.
  REQUIRE(in_degree(1, g) >= 1);
}
