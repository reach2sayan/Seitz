// The generation layer's two algorithmic pieces against brute-force
// references: the lazy Wyckoff-assignment enumeration against the plain
// backtracker it replaced (equal as sets), and the bucketed distance check
// against the all-pairs scan.

#include <cppcrystal/data/element_data.hpp>
#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/generate/distance_check.hpp>
#include <cppcrystal/generate/generator.hpp>
#include <cppcrystal/group/space_group.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <random>
#include <set>
#include <stdexcept>
#include <vector>

using namespace cppcrystal;

namespace {

using cppcrystal::test::must;

// An assignment as a comparable value: sorted (type, position index) pairs.
using Canonical = std::vector<std::pair<int, int>>;

Canonical canonical(generate::Assignment<group::Wyckoff> const &assignment,
                    std::span<group::Wyckoff const> positions) {
  Canonical out;
  for (auto const &p : assignment) {
    out.emplace_back(p.type, static_cast<int>(p.position - positions.data()));
  }
  std::ranges::sort(out);
  return out;
}

// The plain depth-first backtracker the generator replaced, kept as the
// reference: no pruning, explicit undo.
struct Reference {
  std::span<group::Wyckoff const> positions;
  std::vector<std::pair<int, int>> elements;
  std::vector<bool> used_special = std::vector<bool>(positions.size(), false);
  generate::Assignment<group::Wyckoff> placements;
  std::set<Canonical> out;

  void recurse(std::size_t elem, std::size_t pos, int remaining) {
    if (remaining == 0) {
      if (elem + 1 == elements.size()) {
        out.insert(canonical(placements, positions));
      } else {
        recurse(elem + 1, 0, elements[elem + 1].second);
      }
      return;
    }
    if (pos >= positions.size()) {
      return;
    }
    auto const &wp = positions[pos];
    int const mult = wp.multiplicity();
    bool const fixed = wp.degrees_of_freedom() == 0;
    int const max_copies =
        fixed ? (used_special[pos] ? 0 : 1) : remaining / mult;
    for (int copies = 0; copies <= max_copies && copies * mult <= remaining;
         ++copies) {
      for (int k = 0; k < copies; ++k) {
        placements.push_back({elements[elem].first, &wp});
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

std::set<Canonical> reference_assignments(group::SpaceGroup const &sg,
                                          generate::Composition const &comp) {
  Reference ref{sg.wyckoffs(), {}};
  for (auto const &[type, count] : comp) {
    if (count > 0) {
      ref.elements.emplace_back(type, count);
    }
  }
  if (ref.elements.empty()) {
    return {};
  }
  ref.recurse(0, 0, ref.elements.front().second);
  return ref.out;
}

std::set<Canonical> enumerated_assignments(group::SpaceGroup const &sg,
                                           generate::Composition const &comp) {
  std::set<Canonical> out;
  for (auto const &assignment :
       generate::enumerate_assignments(sg.wyckoffs(), comp)) {
    out.insert(canonical(assignment, sg.wyckoffs()));
  }
  return out;
}

bool reference_distances_valid(Cell const &cell,
                               generate::DistanceTolerance tol) {
  auto const radius = [&](Index i) {
    return data::covalent_radius(cell.type(i)).value_or(tol.fallback_radius);
  };
  for (Index i = 0; i < cell.size(); ++i) {
    for (Index j = i; j < cell.size(); ++j) {
      double const d = generate::minimum_image_distance(
          cell.position(i), cell.position(j), cell.lattice().matrix(),
          cell.periodicity(),
          i == j ? generate::Images::nontrivial : generate::Images::all);
      if (d < tol.scale * (radius(i) + radius(j))) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

TEST_CASE("enumerate_assignments equals the backtracking reference as sets",
          "[generate]") {
  struct Case {
    int number;
    generate::Composition comp;
  };
  std::vector<Case> const cases{
      {225, {{11, 4}, {17, 4}}},         // NaCl on Fm-3m
      {221, {{55, 1}, {22, 1}, {8, 3}}}, // perovskite on Pm-3m
      {194, {{30, 2}, {16, 2}}},         // wurtzite-like on P6_3/mmm
      {62, {{26, 4}, {16, 8}}},          // Pnma
      {2, {{6, 6}}},                     // P-1, general position only
      {1, {{6, 3}, {8, 2}}},
      {227, {{14, 8}}},         // diamond on Fd-3m
      {136, {{22, 2}, {8, 4}}}, // rutile
      {225, {{11, 3}}},         // incompatible: 3 on an fcc lattice
      {230, {{6, 1}}},          // incompatible
  };
  for (auto const &c : cases) {
    INFO("space group " << c.number);
    auto const *sg =
        must(group::SpaceGroup::from_number(GroupFamily::space, c.number));
    auto const expected = reference_assignments(*sg, c.comp);
    CHECK(enumerated_assignments(*sg, c.comp) == expected);
    CHECK(generate::Generator{*sg}.compatible(c.comp) == !expected.empty());
  }
}

TEST_CASE("a bounded enumeration stops early", "[generate]") {
  auto const *sg = must(group::SpaceGroup::from_number(GroupFamily::space, 1));
  generate::Composition const comp{{6, 4}, {8, 4}};
  generate::Generator const gen{*sg};
  CHECK(gen.assignments(comp, 1).size() == 1);
  CHECK(gen.assignments(comp, 0).empty());
}

TEST_CASE("distances_valid agrees with the all-pairs scan", "[generate]") {
  std::mt19937 rng(31);
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> wide(-4.0, 4.0);
  constexpr CellPeriodicity layer{AxisKind::periodic, AxisKind::periodic,
                                  AxisKind::aperiodic};
  int agreements = 0;
  for (int trial = 0; trial < 60; ++trial) {
    // A lattice from 3 to 12 angstrom on a side, mildly skewed, and enough
    // atoms of mixed radius that roughly half the trials violate.
    Matrix3d lattice = Matrix3d::Identity() * (3.0 + 9.0 * unit(rng));
    lattice(0, 1) = 2.0 * unit(rng);
    lattice(1, 2) = 1.5 * unit(rng);
    CellPeriodicity const periodicity = trial % 3 == 0   ? all_periodic()
                                        : trial % 3 == 1 ? layer
                                                         : none_periodic();
    Index const n = 4 + static_cast<Index>(rng() % 20);
    Positions pos(n, 3);
    Types types;
    for (Index i = 0; i < n; ++i) {
      for (Index k = 0; k < 3; ++k) {
        pos(i, k) =
            periodicity[static_cast<std::size_t>(k)] == AxisKind::periodic
                ? unit(rng)
                : wide(rng);
      }
      types.push_back(std::array{1, 6, 8, 26, 55}[rng() % 5]);
    }
    Cell const cell(Lattice{lattice}, pos, types, periodicity);
    generate::DistanceTolerance const tol{.scale = 0.4 + 0.6 * unit(rng)};
    bool const expected = reference_distances_valid(cell, tol);
    CHECK(generate::distances_valid(cell, tol) == expected);
    agreements += expected ? 1 : 0;
  }
  // Both outcomes must actually be exercised.
  CHECK(agreements > 5);
  CHECK(agreements < 55);
}
