#pragma once

#include <cppcrystal/data/subgroup_relations.hpp>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace cppcrystal::group {

inline constexpr int kNumSpaceGroups = 230;

// A maximal-subgroup (or minimal-supergroup) relation: the related space-group
// number and the prime index of the relation.
struct SubgroupRelation {
  int number;
  int index;
};

namespace detail {

// The edge list grouped by one endpoint: `edges` sorted by key, with `offsets`
// marking each key's block. A counting sort at compile time, so both adjacency
// directions and the reachability closure below are constant-initialised — the
// runtime does no graph construction at all.
struct Adjacency {
  std::array<SubgroupRelation, data::kTSubgroupRelations.size()> edges{};
  std::array<int, kNumSpaceGroups + 2> offsets{};

  [[nodiscard]] constexpr std::span<SubgroupRelation const>
  operator[](int key) const noexcept {
    if (key < 1 || key > kNumSpaceGroups) {
      return {};
    }
    auto const k = static_cast<std::size_t>(key);
    return {edges.data() + offsets[k], edges.data() + offsets[k + 1]};
  }
};

// `Endpoint` picks which end of the relation groups it (super for out-edges,
// sub for in-edges); `Other` is the number stored alongside the index.
template <int data::TSubgroupRelation::*Endpoint,
          int data::TSubgroupRelation::*Other>
[[nodiscard]] consteval Adjacency group_by() {
  Adjacency out{};
  for (auto const &rel : data::kTSubgroupRelations) {
    ++out.offsets[static_cast<std::size_t>(rel.*Endpoint) + 1];
  }
  for (std::size_t k = 1; k < out.offsets.size(); ++k) {
    out.offsets[k] += out.offsets[k - 1];
  }
  auto cursor = out.offsets;
  for (auto const &rel : data::kTSubgroupRelations) {
    auto const k = static_cast<std::size_t>(rel.*Endpoint);
    out.edges[static_cast<std::size_t>(cursor[k]++)] =
        SubgroupRelation{rel.*Other, rel.index};
  }
  return out;
}

inline constexpr Adjacency kSubgroupsOf =
    group_by<&data::TSubgroupRelation::super, &data::TSubgroupRelation::sub>();
inline constexpr Adjacency kSupergroupsOf =
    group_by<&data::TSubgroupRelation::sub, &data::TSubgroupRelation::super>();

using ReachabilityRow = std::bitset<kNumSpaceGroups + 1>;

// Transitive closure of the maximal-subgroup relation, by Warshall over rows:
// ~53k row-ORs, all at compile time, so is_subgroup is a bit test.
[[nodiscard]] consteval std::array<ReachabilityRow, kNumSpaceGroups + 1>
reachability() {
  std::array<ReachabilityRow, kNumSpaceGroups + 1> reach{};
  for (auto const &rel : data::kTSubgroupRelations) {
    reach[static_cast<std::size_t>(rel.super)].set(
        static_cast<std::size_t>(rel.sub));
  }
  for (std::size_t k = 1; k <= kNumSpaceGroups; ++k) {
    for (std::size_t i = 1; i <= kNumSpaceGroups; ++i) {
      if (reach[i].test(k)) {
        reach[i] |= reach[k];
      }
    }
  }
  return reach;
}

inline constexpr auto kReachable = reachability();

} // namespace detail

// The translationengleiche (lattice-preserving) maximal-subgroup graph of the
// 230 space groups (the t-subgroup part of the group-subgroup relations).
//
// Derived without external maximal-subgroup tables: for each space group the
// maximal subgroups of its point group are enumerated, and the operations that
// survive each one form a t-subgroup, identified back to a space-group number
// by the determination from a set of operations. Klassengleiche
// (cell-multiplying) subgroups are intentionally NOT included. That derivation
// is the offline path (derive_t_subgroup_edges) that produced the baked table;
// the graph itself is constant-initialised from it, so there is no runtime
// construction, no shared state and nothing to warm up.
class SubgroupGraph {
public:
  // Maximal t-subgroups of `number` (out-edges): each a space-group number with
  // the prime index of the relation. Empty for P1 (number 1) and out-of-range.
  [[nodiscard]] static constexpr std::span<SubgroupRelation const>
  maximal_subgroups(int number) noexcept {
    return detail::kSubgroupsOf[number];
  }

  // Minimal t-supergroups of `number` (in-edges): the groups of which `number`
  // is a maximal t-subgroup.
  [[nodiscard]] static constexpr std::span<SubgroupRelation const>
  minimal_supergroups(int number) noexcept {
    return detail::kSupergroupsOf[number];
  }

  // True iff `sub` is a translationengleiche subgroup of `super` — reachable
  // from `super` by a chain of maximal t-subgroup steps. True when sub ==
  // super. A bit test in the compile-time closure.
  [[nodiscard]] static constexpr bool is_subgroup(int sub, int super) noexcept {
    if (super < 1 || super > kNumSpaceGroups || sub < 1 ||
        sub > kNumSpaceGroups) {
      return false;
    }
    return sub == super ||
           detail::kReachable[static_cast<std::size_t>(super)].test(
               static_cast<std::size_t>(sub));
  }

  // A chain of space-group numbers descending from `super` to `sub` along
  // maximal t-subgroup steps (super first, sub last); std::nullopt if `sub` is
  // not a t-subgroup of `super`. A symmetry-breaking transition path.
  [[nodiscard]] static std::optional<std::vector<int>> path(int super, int sub);
};

// One derived edge: `sub` is a maximal translationengleiche subgroup of `super`
// at the given prime index.
struct TSubgroupEdge {
  int super;
  int sub;
  int index;
};

// Re-derive the whole t-subgroup edge list from scratch (enumerate each space
// group's maximal point-subgroups and identify the surviving operations via the
// determination). This is the SLOW offline path used only to regenerate the
// baked table data/subgroup_relations.hpp (see
// tools/generate_subgroup_relations.cpp).
[[nodiscard]] std::vector<TSubgroupEdge> derive_t_subgroup_edges();

} // namespace cppcrystal::group
