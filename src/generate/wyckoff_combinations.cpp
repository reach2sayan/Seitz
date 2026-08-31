#include <cppcrystal/generate/wyckoff_combinations.hpp>

namespace cppcrystal::generate {

std::vector<WyckoffCombination>
list_wyckoff_combinations(group::SpaceGroup const &sg, Composition const &comp,
                          std::size_t max_combinations) {
  return enumerate_assignments(sg.wyckoffs(), comp, max_combinations);
}

bool check_compatible(group::SpaceGroup const &sg, Composition const &comp) {
  return !list_wyckoff_combinations(sg, comp, 1).empty();
}

} // namespace cppcrystal::generate
