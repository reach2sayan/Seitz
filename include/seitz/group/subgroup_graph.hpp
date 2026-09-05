#pragma once

#include <seitz/data/detail/lookup.hpp>
#include <seitz/data/subgroup_relations.hpp>

#include <array>
#include <bitset>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::group {

inline constexpr int kNumSpaceGroups = 230;

// A maximal-subgroup (or minimal-supergroup) relation: the related space-group
// number and the prime index of the relation.
struct SubgroupRelation {
  int number;
  int index;
};

namespace detail {

// The edge list grouped by one endpoint, via the consteval counting sort the
// data catalogs use (data/detail/lookup.hpp). Both directions and the closure
// below are constant-initialised: no runtime graph construction.
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

// Transitive closure of the maximal-subgroup relation, Warshall over rows:
// ~53k row-ORs at compile time, so is_subgroup is a bit test.
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
// 230 space groups.
//
// Derived without external tables: per group, the maximal subgroups of its
// point group are enumerated and the operations surviving each -- a t-subgroup
// -- identified back to a number by determination. Klassengleiche
// (cell-multiplying) subgroups are deliberately excluded. That derivation
// (derive_t_subgroup_edges) is offline; the table it baked is
// constant-initialised here, so nothing is built or warmed up at runtime.
class SubgroupGraph {
public:
  // Maximal t-subgroups of `number` (out-edges), each with the prime index of
  // the relation. Empty for P1 and out-of-range.
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

  // `sub` reachable from `super` by a chain of maximal t-subgroup steps; true
  // when sub == super. A bit test in the compile-time closure.
  [[nodiscard]] static constexpr bool is_subgroup(int sub, int super) noexcept {
    if (super < 1 || super > kNumSpaceGroups || sub < 1 ||
        sub > kNumSpaceGroups) {
      return false;
    }
    return sub == super ||
           detail::kReachable[static_cast<std::size_t>(super)].test(
               static_cast<std::size_t>(sub));
  }

  // A chain super -> ... -> sub along maximal t-subgroup steps (super first),
  // nullopt if `sub` is no t-subgroup of `super`: a symmetry-breaking path.
  [[nodiscard]] static std::optional<std::vector<int>> path(int super, int sub);
};

// One derived edge: `sub` is a maximal translationengleiche subgroup of `super`
// at the given prime index.
struct TSubgroupEdge {
  int super;
  int sub;
  int index;
};

// Re-derive the whole edge list (maximal point-subgroups + determination of
// the survivors). SLOW; offline only, to regenerate tools/subgroup_relations.csv
// (tools/generate_subgroup_relations.cpp), which the build transcribes.
[[nodiscard]] std::vector<TSubgroupEdge> derive_t_subgroup_edges();

} // namespace seitz::group

#pragma GCC visibility pop
