#pragma once

#include <boost/container/flat_map.hpp>

#include <seitz/core/mdspan.hpp>
#include <seitz/generate/distance_check.hpp>

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

namespace seitz::generate {

// Atom type -> atom count.
using Composition = boost::container::flat_map<int, int>;

// The interface a Wyckoff position must offer to be assignable.
template <class W>
concept WyckoffLike = requires(W const &w) {
  { w.multiplicity() } -> std::same_as<int>;
  { w.degrees_of_freedom() } -> std::same_as<int>;
};

// An atom type placed on a chosen Wyckoff position (non-owning pointer into
// the group, which must outlive it). `coordinate` is the fractional generating
// point of the orbit: set by the caller to pin it, filled in by the generator
// on the assignment it returns.
template <WyckoffLike W> struct Placed {
  int type;
  W const *position{};
  std::optional<Vector3d> coordinate = std::nullopt;
};

// One complete Wyckoff assignment of a composition.
template <WyckoffLike W> using Assignment = std::vector<Placed<W>>;

// Which Wyckoff positions a structure may be built on. `general_only` keeps
// generic coordinates, so no special site adds accidental symmetry (a layer
// group's atoms pinned in one plane gain a horizontal mirror); it requires
// N_atoms to be a multiple of the general multiplicity.
enum class Placement { any, general_only };

// An atom pre-assigned to a Wyckoff position by letter, and optionally to a
// generating coordinate on it (projected onto the position's locus).
struct FixedSite {
  int type;
  char letter;
  std::optional<Vector3d> coordinate = std::nullopt;
};

// What a Generator may vary while searching for a structure.
struct GenerateOptions {
  // Multiplies the element-aware size estimate.
  double scale = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / lattice resampling attempts per Wyckoff assignment.
  int attempts_per_combination = 50;
  DistanceTolerance distance = {}; // min-distance acceptance limit
  Placement placement = Placement::any;
  // Generate into this lattice instead of a random one. Its metric must be
  // invariant under the group's operations (e_incompatible_lattice).
  std::optional<Lattice> lattice = std::nullopt;
  // Atoms placed before the search; counted against the composition.
  std::vector<FixedSite> sites = {};
};

namespace detail {

// Bound on positions per group (27 in 3D, fewer for layer and rod).
constexpr std::size_t kMaxPositions = 64;
using UsedSpecial = std::bitset<kMaxPositions>;

// Fixed data of one enumeration: the elements to place, plus
// reachable_[p, r] = can positions p.. supply exactly r atoms. The table
// ignores the once-only rule on fixed positions, so it over-approximates --
// a valid prune.
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
    // Named, not inline as `table[...size(), 0]`: MSVC applies its
    // discarded-[[nodiscard]] check to non-final operands of a
    // multidimensional subscript, so span::size() there raises C4834.
    auto const last = static_cast<Index>(positions.size());
    auto const rows = positions.size() + 1;
    auto const cols = static_cast<std::size_t>(ctx.max_count) + 1;
    ctx.reachable_.assign(rows * cols, 0);
    md::matrix_view<std::uint8_t> table(ctx.reachable_.data(),
                                        static_cast<Index>(rows),
                                        static_cast<Index>(cols));
    table[last, 0] = 1;
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

// Depth-first: per element, choose each position's copy count (0 or 1 with no
// DOF, any number with free coordinates) so the multiplicities sum to that
// element's count. `used_special` enforces one use per no-DOF position -- it is
// a single fixed orbit. Yields the shared `placements` buffer.
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

// Every valid Wyckoff assignment of `comp` on `positions`, lazily, each
// beginning with the `fixed` placements (which point into `positions` and are
// deducted from `comp`). The yielded reference points at a reused buffer: copy
// what you keep. Bound it with views::take.
template <WyckoffLike W>
[[nodiscard]] std::generator<Assignment<W> const &>
enumerate_assignments(std::span<W const> positions, Composition comp,
                      Assignment<W> fixed = {}) {
  if (positions.size() > detail::kMaxPositions) {
    co_return;
  }
  detail::UsedSpecial used;
  for (auto const &placed : fixed) {
    auto const index = static_cast<std::size_t>(placed.position - positions.data());
    if (index >= positions.size()) {
      co_return; // not a position of this group
    }
    comp[placed.type] -= placed.position->multiplicity();
    if (placed.position->degrees_of_freedom() == 0) {
      if (used[index]) {
        co_return; // one fixed orbit, placed twice
      }
      used.set(index);
    }
  }
  if (std::ranges::any_of(comp | std::views::values,
                          [](int count) { return count < 0; })) {
    co_return; // the fixed sites overshoot the composition
  }
  auto const ctx = detail::AssignmentContext<W>::of(positions, comp);
  if (ctx.elements.empty()) {
    if (!fixed.empty()) {
      co_yield fixed; // the fixed sites are the whole structure
    }
    co_return;
  }
  Assignment<W> placements = std::move(fixed);
  co_yield std::ranges::elements_of(detail::walk(
      ctx, placements, used, 0, 0, ctx.elements.front().second));
}

// Whether `comp` has at least one assignment on `positions`.
template <WyckoffLike W>
[[nodiscard]] bool assignable(std::span<W const> positions,
                              Composition const &comp) {
  auto assignments = enumerate_assignments(positions, comp);
  return assignments.begin() != assignments.end();
}

} // namespace seitz::generate

#pragma GCC visibility pop
