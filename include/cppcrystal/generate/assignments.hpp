#pragma once

#include <boost/container/flat_map.hpp>

#include <cppcrystal/core/mdspan.hpp>
#include <cppcrystal/generate/distance_check.hpp>

#include <algorithm>
#include <bitset>
#include <concepts>
#include <cstdint>
#include <generator>
#include <iterator>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

// The shared vocabulary of random-structure generation (3D, layer, rod,
// cluster): Wyckoff assignments of a composition and their enumeration. The
// search over them is generate::Generator.

#pragma GCC visibility push(default)

namespace cppcrystal::generate {

// Atom type -> atom count.
using Composition = boost::container::flat_map<int, int>;

// The interface a Wyckoff position must offer to be assignable.
template <class W>
concept WyckoffLike = requires(W const &w) {
  { w.multiplicity() } -> std::same_as<int>;
  { w.degrees_of_freedom() } -> std::same_as<int>;
};

// An atom type placed on a chosen Wyckoff position (non-owning pointer into
// the group, which must outlive it).
template <WyckoffLike W> struct Placed {
  int type;
  W const *position{};
};

// One complete Wyckoff assignment of a composition.
template <WyckoffLike W> using Assignment = std::vector<Placed<W>>;

// Which Wyckoff positions a structure may be built on. `general_only` keeps
// generic coordinates only, so the exact target group is realized without the
// accidental extra symmetry a special (fixed/high-symmetry) site can introduce
// — notably for layer groups, where atoms confined to a special site in one
// plane gain a horizontal mirror. It requires the total atom count to be a
// multiple of the general-position multiplicity.
enum class Placement { any, general_only };

// What a Generator may vary while searching for a structure.
struct GenerateOptions {
  // Multiplies the element-aware size estimate (cell volume, cluster metric,
  // or rod repeat length — whichever the family uses).
  double scale = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / lattice resampling attempts per Wyckoff assignment.
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion for a generated structure.
  DistanceTolerance distance = {};
  Placement placement = Placement::any;
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
    AssignmentContext ctx{
        .positions = positions,
    };
    ctx.elements = {std::from_range,
                    comp | std::views::filter([](auto const &entry) {
                      return entry.second > 0;
                    })};
    if (ctx.elements.empty()) {
      return ctx;
    }
    ctx.max_count = std::ranges::max(ctx.elements | std::views::values);
    // Named rather than written inline as `table[...size(), 0]`: MSVC applies
    // its discarded-[[nodiscard]] check to the non-final operands of a
    // multidimensional subscript, so a call there -- span::size() is
    // [[nodiscard]] in that standard library -- raises a spurious C4834.
    auto const last = static_cast<Index>(positions.size());
    auto const rows = positions.size() + 1;
    auto const cols = static_cast<std::size_t>(ctx.max_count) + 1;
    ctx.reachable_.assign(rows * cols, 0);
    md::matrix_view<std::uint8_t> table(ctx.reachable_.data(),
                                        static_cast<Index>(rows),
                                        static_cast<Index>(cols));
    table[last, 0] = 1;
    // Backwards over the positions: take zero copies of p, or one copy and
    // stay on p (free coordinates) / move past it (fixed position).
    for (auto p = last; p-- > 0;) {
      auto const &wp = positions[static_cast<std::size_t>(p)];
      auto const mult = static_cast<Index>(wp.multiplicity());
      Index const after_one = wp.degrees_of_freedom() == 0 ? p + 1 : p;
      for (Index r = 0; r <= ctx.max_count; ++r) {
        table[p, r] =
            table[p + 1, r] || (r >= mult && table[after_one, r - mult]);
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
      co_yield std::ranges::elements_of(walk(ctx, placements, used_special,
                                             elem + 1, 0,
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
  int const max_copies = fixed ? (used_special[pos] ? 0 : 1) : remaining / mult;

  for (int copies = 0; copies <= max_copies && copies * mult <= remaining;
       ++copies) {
    placements.insert(placements.end(), static_cast<std::size_t>(copies),
                      Placed<W>{ctx.elements[elem].first, &wp});
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
  co_yield std::ranges::elements_of(
      detail::walk(ctx, placements, {}, 0, 0, ctx.elements.front().second));
}

// Whether `comp` has at least one assignment on `positions`.
template <WyckoffLike W>
[[nodiscard]] bool assignable(std::span<W const> positions,
                              Composition const &comp) {
  auto assignments = enumerate_assignments(positions, comp);
  return assignments.begin() != assignments.end();
}

} // namespace cppcrystal::generate

#pragma GCC visibility pop
