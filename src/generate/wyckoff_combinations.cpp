#include <cppcrystal/generate/wyckoff_combinations.hpp>

#include <span>
#include <vector>

namespace cppcrystal::generate {

namespace {

// One element of the composition.
struct Element {
  int type;
  int count;
};

// Depth-first enumerator of Wyckoff assignments. It walks the elements in order;
// for each element it chooses how many copies of each Wyckoff position to use
// (0 or 1 for a fixed/no-DOF position, any number for a position with free
// coordinates), so that the chosen multiplicities sum to the element's count.
// A `used_special` flag per position enforces the "a no-DOF position is one
// fixed orbit, usable at most once across the whole structure" rule.
struct Enumerator {
  std::span<group::WyckoffPosition const> positions;
  std::vector<Element> const &elements;
  std::size_t max_combinations;

  std::vector<bool> used_special;
  std::vector<WyckoffCombination::Placement> placements;
  std::vector<WyckoffCombination> out;

  [[nodiscard]] bool full() const { return out.size() >= max_combinations; }

  // Choose copies of positions[pos] onward to make `remaining` atoms of the
  // current element, then move on to the next element.
  void recurse(std::size_t elem, std::size_t pos, int remaining) {
    if (full()) {
      return;
    }
    if (remaining == 0) {
      if (elem + 1 == elements.size()) {
        out.push_back(WyckoffCombination{placements});
      } else {
        recurse(elem + 1, 0, elements[elem + 1].count);
      }
      return;
    }
    if (pos >= positions.size()) {
      return; // ran out of positions before filling the element
    }

    group::WyckoffPosition const &wp = positions[pos];
    int const mult = wp.multiplicity();
    bool const fixed = wp.degrees_of_freedom() == 0;
    int const max_copies =
        fixed ? (used_special[pos] ? 0 : 1) : remaining / mult;

    for (int copies = 0; copies <= max_copies && copies * mult <= remaining;
         ++copies) {
      if (full()) {
        return;
      }
      for (int k = 0; k < copies; ++k) {
        placements.push_back({elements[elem].type, &wp});
      }
      if (fixed && copies == 1) {
        used_special[pos] = true;
      }

      recurse(elem, pos + 1, remaining - copies * mult);

      if (fixed && copies == 1) {
        used_special[pos] = false;
      }
      for (int k = 0; k < copies; ++k) {
        placements.pop_back();
      }
    }
  }
};

} // namespace

std::vector<WyckoffCombination>
list_wyckoff_combinations(group::SpaceGroup const &sg, Composition const &comp,
                          std::size_t max_combinations) {
  std::vector<Element> elements;
  elements.reserve(comp.size());
  for (auto const &[type, count] : comp) {
    if (count > 0) {
      elements.push_back({type, count});
    }
  }
  if (elements.empty() || max_combinations == 0) {
    return {};
  }

  Enumerator e{sg.wyckoffs(),
               elements,
               max_combinations,
               std::vector<bool>(sg.wyckoffs().size(), false),
               {},
               {}};
  e.recurse(0, 0, elements.front().count);
  return std::move(e.out);
}

bool check_compatible(group::SpaceGroup const &sg, Composition const &comp) {
  return !list_wyckoff_combinations(sg, comp, 1).empty();
}

} // namespace cppcrystal::generate
