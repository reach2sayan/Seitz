#pragma once

#include <memory>
#include <optional>
#include <vector>

namespace cppcrystal::group {

// A maximal-subgroup (or minimal-supergroup) relation: the related space-group
// number and the prime index of the relation.
struct SubgroupRelation {
  int number;
  int index;
};

// The translationengleiche (lattice-preserving) maximal-subgroup graph of the
// 230 space groups (the t-subgroup part of the group-subgroup relations).
//
// Derived without external maximal-subgroup tables. For each space group the
// maximal subgroups of its point group are enumerated; the
// operations that survive each one form a t-subgroup, which is identified back
// to a space-group number by the (oracle-validated) determination from a set of
// operations (spacegroup_type_from_symmetry). Klassengleiche (cell-multiplying)
// subgroups are intentionally NOT included.
//
// The relation lattice is held in a Boost.Graph and built once, lazily, on
// first access. Construction runs a determination per candidate subgroup, so
// the first call to instance() does real work; subsequent queries are cheap.
class SubgroupGraph {
public:
  SubgroupGraph(SubgroupGraph const &) = delete;
  SubgroupGraph &operator=(SubgroupGraph const &) = delete;
  ~SubgroupGraph();

  // The shared, lazily-built graph.
  [[nodiscard]] static SubgroupGraph const &instance();

  // Maximal t-subgroups of `number` (out-edges): each a space-group number with
  // the prime index of the relation. Empty for P1 (number 1) and out-of-range.
  [[nodiscard]] std::vector<SubgroupRelation>
  maximal_subgroups(int number) const;

  // Minimal t-supergroups of `number` (in-edges): the groups of which `number`
  // is a maximal t-subgroup.
  [[nodiscard]] std::vector<SubgroupRelation>
  minimal_supergroups(int number) const;

  // True iff `sub` is a translationengleiche subgroup of `super` — reachable
  // from `super` by a chain of maximal t-subgroup steps. True when sub ==
  // super.
  [[nodiscard]] bool is_subgroup(int sub, int super) const;

  // A chain of space-group numbers descending from `super` to `sub` along
  // maximal t-subgroup steps (super first, sub last); std::nullopt if `sub` is
  // not a t-subgroup of `super`. A symmetry-breaking transition path.
  [[nodiscard]] std::optional<std::vector<int>> path(int super, int sub) const;

private:
  SubgroupGraph();

  struct Impl;
  std::unique_ptr<Impl> impl_;
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
// tools/generate_subgroup_relations .cpp); the runtime SubgroupGraph loads that
// table instead of calling this.
[[nodiscard]] std::vector<TSubgroupEdge> derive_t_subgroup_edges();

} // namespace cppcrystal::group
