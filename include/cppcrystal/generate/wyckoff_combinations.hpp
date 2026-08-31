#pragma once

#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/group/space_group.hpp>
#include <cppcrystal/group/wyckoff_position.hpp>

#include <vector>

namespace cppcrystal::generate {

// One complete Wyckoff assignment on a space (or layer) group.
using WyckoffCombination = Assignment<group::WyckoffPosition>;

// Whether `comp` can be placed on `sg` at all.
[[nodiscard]] bool check_compatible(group::SpaceGroup const &sg,
                                    Composition const &comp);

// All valid Wyckoff assignments for `comp` on `sg` (see
// generate::enumerate_assignments for the rules and the cap semantics).
[[nodiscard]] std::vector<WyckoffCombination>
list_wyckoff_combinations(group::SpaceGroup const &sg, Composition const &comp,
                          std::size_t max_combinations = 1000);

} // namespace cppcrystal::generate
