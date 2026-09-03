#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/mdspan.hpp>
#include <cppcrystal/generate/distance_check.hpp>

#include <algorithm>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <format>
#include <generator>
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
template <typename Realize, WyckoffLike W>
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

// A structure never has more Wyckoff positions than this (27 in 3D, fewer
// for layer and rod groups); the once-only bookkeeping of fixed positions is a
// bitset of that width, passed by value down the search.
constexpr std::size_t kMaxPositions = 64;
using UsedSpecial = std::bitset<kMaxPositions>;

// The fixed data of one enumeration: the elements to place and, per (position,
// remainder), whether the positions from there on can supply exactly that many
// atoms. The table ignores the once-only rule on fixed positions, so it
// over-approximates and is therefore a valid prune of the search.
template <WyckoffLike W> struct AssignmentContext {
  std::span<W const> positions;
  std::vector<std::pair<int, int>> elements{}; // (type, count), count > 0
  std::vector<std::uint8_t> reachable_{};      // (positions + 1) x (max + 1)
  int max_count = 0;

  [[nodiscard]] static AssignmentContext of(std::span<W const> positions,
                                            Composition const &comp) {
    AssignmentContext ctx{.positions = positions,};
    std::ranges::copy(comp | std::views::filter([](auto const &entry) {
                        return entry.second > 0;
                      }),
                      std::back_inserter(ctx.elements));
    if (ctx.elements.empty()) {
      return ctx;
    }
    ctx.max_count = std::ranges::max(ctx.elements | std::views::values);
    auto const rows = positions.size() + 1;
    auto const cols = static_cast<std::size_t>(ctx.max_count) + 1;
    ctx.reachable_.assign(rows * cols, 0);
    md::matrix_view<std::uint8_t> table(ctx.reachable_.data(),
                                        static_cast<Index>(rows),
                                        static_cast<Index>(cols));
    table[static_cast<Index>(positions.size()), 0] = 1;
    // Backwards over the positions: take zero copies of p, or one copy and
    // stay on p (free coordinates) / move past it (fixed position).
    for (auto p = static_cast<Index>(positions.size()); p-- > 0;) {
      auto const &wp = positions[static_cast<std::size_t>(p)];
      auto const mult = static_cast<Index>(wp.multiplicity());
      Index const after_one = wp.degrees_of_freedom() == 0 ? p + 1 : p;
      for (Index r = 0; r <= ctx.max_count; ++r) {
        table[p, r] = table[p + 1, r] ||
                      (r >= mult && table[after_one, r - mult]);
      }
    }
    return ctx;
  }

  [[nodiscard]] bool reachable(std::size_t pos, int remaining) const noexcept {
    md::matrix_view<std::uint8_t const> const table(
        reachable_.data(), static_cast<Index>(positions.size() + 1),
        static_cast<Index>(max_count) + 1);
    return table[static_cast<Index>(pos), remaining] != 0;
  }
};

// Depth-first walk of the assignments. Element by element, choose how many
// copies of each position to use (0 or 1 for a fixed/no-DOF position, any
// number for a position with free coordinates) so the chosen multiplicities
// sum to the element's count. `used_special` enforces the "a no-DOF position
// is one fixed orbit, usable at most once across the whole structure" rule.
// Yields a reference to the shared `placements` buffer for each complete
// assignment.
template <WyckoffLike W>
std::generator<Assignment<W> const &>
walk(AssignmentContext<W> const &ctx, Assignment<W> &placements,
     UsedSpecial used_special, std::size_t elem, std::size_t pos,
     int remaining) {
  if (remaining == 0) {
    if (elem + 1 == ctx.elements.size()) {
      co_yield placements;
    } else {
      co_yield std::ranges::elements_of(
          walk(ctx, placements, used_special, elem + 1, 0,
               ctx.elements[elem + 1].second));
    }
    co_return;
  }
  if (pos >= ctx.positions.size() || !ctx.reachable(pos, remaining)) {
    co_return; // the remaining positions cannot fill the element
  }

  W const &wp = ctx.positions[pos];
  int const mult = wp.multiplicity();
  bool const fixed = wp.degrees_of_freedom() == 0;
  int const max_copies =
      fixed ? (used_special[pos] ? 0 : 1) : remaining / mult;

  for (int copies = 0; copies <= max_copies && copies * mult <= remaining;
       ++copies) {
    placements.insert(placements.end(), static_cast<std::size_t>(copies),
                      Placement<W>{ctx.elements[elem].first, &wp});
    UsedSpecial used = used_special;
    if (fixed && copies == 1) {
      used.set(pos);
    }
    co_yield std::ranges::elements_of(
        walk(ctx, placements, used, elem, pos + 1, remaining - copies * mult));
    placements.erase(placements.end() - copies, placements.end());
  }
}

} // namespace detail

// Every valid Wyckoff assignment of `comp` on `positions`, lazily and in
// depth-first order. Each yielded reference points at a buffer reused for
// the next assignment: copy what you keep. Compose with views::take for a
// bounded enumeration; a composition with no assignment yields nothing
// (usually without walking the tree, thanks to the reachability prune).
template <WyckoffLike W>
[[nodiscard]] std::generator<Assignment<W> const &>
enumerate_assignments(std::span<W const> positions, Composition comp) {
  auto const ctx = detail::AssignmentContext<W>::of(positions, comp);
  if (ctx.elements.empty() || positions.size() > detail::kMaxPositions) {
    co_return;
  }
  Assignment<W> placements;
  co_yield std::ranges::elements_of(detail::walk(
      ctx, placements, {}, 0, 0, ctx.elements.front().second));
}

// Whether `comp` has at least one assignment on `positions`.
template <WyckoffLike W>
[[nodiscard]] bool assignable(std::span<W const> positions,
                              Composition const &comp) {
  auto assignments = enumerate_assignments(positions, comp);
  return assignments.begin() != assignments.end();
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

// Assignments considered per search: enough for any real composition, a
// bound for pathological ones.
constexpr std::size_t kMaxAssignments = 1000;

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
  auto assignments = std::ranges::to<std::vector<Assignment<W>>>(
      enumerate_assignments(positions, comp) |
      std::views::take(kMaxAssignments));
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
