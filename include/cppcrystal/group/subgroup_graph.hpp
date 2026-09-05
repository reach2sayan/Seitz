#pragma once

#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/subgroup_relations.hpp>

#include <array>
#include <bitset>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#pragma GCC visibility push(default)

namespace cppcrystal::group {

inline constexpr int kNumSpaceGroups = 230;

// A maximal-subgroup (or minimal-supergroup) relation: the related space-group
// number and the prime index of the relation.
struct SubgroupRelation {
  int number;
  int index;
};

namespace detail {

// The edge list grouped by one endpoint: the same consteval counting sort the
// data catalogs are indexed with (data/detail/lookup.hpp), over the baked
// relation table. Both adjacency directions and the reachability closure below
// are constant-initialised, so the runtime does no graph construction at all.
inline constexpr std::size_t kNumRelations = data::kTSubgroupRelations.size();
inline constexpr std::size_t kNumKeys =
    static_cast<std::size_t>(kNumSpaceGroups) + 1;

using Adjacency =
    data::detail::BucketIndex<SubgroupRelation, kNumRelations, kNumKeys>;

[[nodiscard]] constexpr data::TSubgroupRelation const &
relation(int id) noexcept {
  return data::kTSubgroupRelations[static_cast<std::size_t>(id) - 1];
}

// `Endpoint` picks which end of the relation groups it (super for out-edges,
// sub for in-edges); `Other` is the number stored alongside the index.
template <int data::TSubgroupRelation::*Endpoint,
          int data::TSubgroupRelation::*Other>
[[nodiscard]] consteval Adjacency group_by() {
  return data::detail::bucket_index<kNumRelations, kNumKeys, SubgroupRelation>(
      [](int id) { return relation(id).*Endpoint; },
      [](int id) {
        return SubgroupRelation{.number = relation(id).*Other,
                                .index = relation(id).index};
      });
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

#pragma GCC visibility pop
