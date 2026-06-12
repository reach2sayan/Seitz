#pragma once

#include <spglib/group/space_group.hpp>
#include <spglib/group/wyckoff_position.hpp>

#include <map>
#include <vector>

namespace spglib::generate {

using Composition = std::map<int, int>;
struct WyckoffCombination {
  struct Placement {
    int type;                                 // atom type
    group::WyckoffPosition const *position{}; // chosen Wyckoff position
  };
  std::vector<Placement> placements;
};

// Whether `comp` can be placed on `sg` at all (PyXtal Group.check_compatible).
[[nodiscard]] bool check_compatible(group::SpaceGroup const &sg,
                                    Composition const &comp);

// All valid Wyckoff assignments for `comp` on `sg` (PyXtal
// Group.list_wyckoff_combinations). A position with no degrees of freedom is a
// single fixed orbit, so it can be used at most once across the whole structure;
// positions with at least one free coordinate may be reused. `max_combinations`
// caps the result to keep enumeration bounded; when the cap is hit, fewer than
// the full set are returned (the caller can detect truncation by size ==
// max_combinations).
[[nodiscard]] std::vector<WyckoffCombination>
list_wyckoff_combinations(group::SpaceGroup const &sg, Composition const &comp,
                          std::size_t max_combinations = 1000);

} // namespace spglib::generate
