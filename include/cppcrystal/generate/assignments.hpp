#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/generate/distance_check.hpp>

#include <algorithm>
#include <concepts>
#include <cstdint>
#include <format>
#include <iterator>
#include <map>
#include <optional>
#include <random>
#include <ranges>
#include <span>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

// The shared vocabulary of the random-structure generators (3D, layer, rod,
// cluster): Wyckoff assignments of a composition, their enumeration, and the
// shuffled assignment/attempt search driver.
namespace cppcrystal::generate {

// Atom type -> atom count.
using Composition = std::map<int, int>;

// The interface a Wyckoff position must offer to be assignable.
template <class W>
concept WyckoffLike = requires(W const &w) {
  { w.multiplicity() } -> std::same_as<int>;
  { w.degrees_of_freedom() } -> std::same_as<int>;
};

// An atom type placed on a chosen Wyckoff position (non-owning pointer into
// the group, which must outlive it).
template <WyckoffLike W> struct Placement {
  int type;
  W const *position{};
};

// One complete Wyckoff assignment of a composition.
template <WyckoffLike W> using Assignment = std::vector<Placement<W>>;

// The candidate type a `realize` callback yields: it returns
// std::optional<R>, and this is that R.
template <class Realize, WyckoffLike W>
using Realized = std::invoke_result_t<Realize &, Assignment<W> const &, int,
                                      std::mt19937_64 &>::value_type;

// Options shared by every random-structure generator.
struct GenerateOptions {
  // Multiplies the element-aware size estimate (cell volume, cluster metric,
  // or rod repeat length — per generator).
  double scale = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / lattice resampling attempts per Wyckoff assignment.
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion for a generated structure.
  DistanceTolerance distance = {};
  // Restrict placement to the general position (generic coordinates only).
  // Useful when the exact target group must be realized without the accidental
  // extra symmetry that special (fixed/high-symmetry) sites can introduce —
  // notably for layer groups, where atoms confined to a special site in one
  // plane gain a horizontal mirror. Requires the total atom count to be a
  // multiple of the general-position multiplicity.
  bool general_position_only = false;
};

namespace detail {

// Depth-first enumerator of Wyckoff assignments. It walks the elements in
// order; for each element it chooses how many copies of each position to use
// (0 or 1 for a fixed/no-DOF position, any number for a position with free
// coordinates), so the chosen multiplicities sum to the element's count.
// `used_special` enforces the "a no-DOF position is one fixed orbit, usable at
// most once across the whole structure" rule.
template <WyckoffLike W> struct AssignmentEnumerator {
  std::span<W const> positions;
  std::vector<std::pair<int, int>> elements; // (type, count)
  std::size_t max_combinations;

  std::vector<bool> used_special = std::vector<bool>(positions.size(), false);
  Assignment<W> placements{};
  std::vector<Assignment<W>> out{};

  [[nodiscard]] bool full() const { return out.size() >= max_combinations; }

  // Choose copies of positions[pos] onward to make `remaining` atoms of the
  // current element, then move on to the next element.
  void recurse(std::size_t elem, std::size_t pos, int remaining) {
    if (full()) {
      return;
    }
    if (remaining == 0) {
      if (elem + 1 == elements.size()) {
        out.push_back(placements);
      } else {
        recurse(elem + 1, 0, elements[elem + 1].second);
      }
      return;
    }
    if (pos >= positions.size()) {
      return; // ran out of positions before filling the element
    }

    W const &wp = positions[pos];
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

} // namespace detail

// All valid Wyckoff assignments for `comp` on `positions`. `max_combinations`
// caps the result to keep enumeration bounded; when the cap is hit, fewer than
// the full set are returned (detectable by size == max_combinations).
template <WyckoffLike W>
[[nodiscard]] std::vector<Assignment<W>>
enumerate_assignments(std::span<W const> positions, Composition const &comp,
                      std::size_t max_combinations = 1000) {
  std::vector<std::pair<int, int>> elements;
  elements.reserve(comp.size());
  std::ranges::copy(comp | std::views::filter([](auto const &entry) {
                      return entry.second > 0;
                    }),
                    std::back_inserter(elements));
  if (elements.empty() || max_combinations == 0) {
    return {};
  }
  detail::AssignmentEnumerator<W> e{positions, std::move(elements),
                                    max_combinations};
  e.recurse(0, 0, e.elements.front().second);
  return std::move(e.out);
}

// Drop every assignment that uses anything but the general position (the last
// position of the group), for general-position-only generation.
template <WyckoffLike W>
void restrict_to_general_position(std::vector<Assignment<W>> &assignments,
                                  std::span<W const> positions) {
  W const *general = &positions.back();
  std::erase_if(assignments, [&](Assignment<W> const &assignment) {
    return std::ranges::any_of(assignment, [&](Placement<W> const &p) {
      return p.position != general;
    });
  });
}

namespace detail {

// The generator pipeline shared by every random-structure entry point:
// enumerate -> optionally restrict to the general position -> shuffle -> for
// each assignment and attempt, ask `realize` for a candidate. `realize` is
// (Assignment<W> const&, int attempt, std::mt19937_64&) -> std::optional<R>;
// the first engaged result wins. `who` names the entry point in errors and
// `group_kind` the group family ("space group", "rod group", ...).
template <WyckoffLike W, class Realize>
  requires std::invocable<Realize &, Assignment<W> const &, int,
                          std::mt19937_64 &>
[[nodiscard]] auto
search_assignments(std::span<W const> positions, Composition const &comp,
                   GenerateOptions const &options, std::string_view who,
                   std::string_view group_kind, Realize &&realize)
    -> Result<Realized<Realize, W>> {
  auto assignments = enumerate_assignments(positions, comp);
  if (assignments.empty()) {
    return leaf::new_error(e_message{std::format(
        "{}: composition is not compatible with the Wyckoff positions of the "
        "requested {}",
        who, group_kind)});
  }

  if (options.general_position_only) {
    restrict_to_general_position(assignments, positions);
    if (assignments.empty()) {
      return leaf::new_error(e_message{std::format(
          "{}: general_position_only requires the atom count to be a multiple "
          "of the general-position multiplicity",
          who)});
    }
  }

  std::mt19937_64 rng(options.seed.value_or(0));
  // Strategic sampling: try the compatible assignments in a random order
  // rather than always the first/one blindly-indexed combination.
  std::ranges::shuffle(assignments, rng);

  for (auto const &assignment : assignments) {
    for (int attempt = 0; attempt < options.attempts_per_combination;
         ++attempt) {
      if (auto candidate = realize(assignment, attempt, rng)) {
        return std::move(*candidate);
      }
    }
  }

  return leaf::new_error(e_message{std::format(
      "{}: no distance-valid structure found within the attempt budget for "
      "any compatible Wyckoff assignment",
      who)});
}

} // namespace detail

} // namespace cppcrystal::generate
