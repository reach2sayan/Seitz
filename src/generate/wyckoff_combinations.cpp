#include <cppcrystal/generate/wyckoff_combinations.hpp>

#include <ranges>

namespace cppcrystal::generate {

std::vector<WyckoffCombination>
list_wyckoff_combinations(group::SpaceGroup const &sg, Composition const &comp,
                          std::size_t max_combinations) {
  return std::ranges::to<std::vector<WyckoffCombination>>(
      enumerate_assignments(sg.wyckoffs(), comp) |
      std::views::take(max_combinations));
}

bool check_compatible(group::SpaceGroup const &sg, Composition const &comp) {
  return assignable(sg.wyckoffs(), comp);
}

} // namespace cppcrystal::generate
