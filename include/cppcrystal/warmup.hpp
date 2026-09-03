#pragma once

#include <future>
#include <utility>

namespace cppcrystal {

// Which per-setting caches to build ahead of first use. group::SpaceGroup::of
// is a flyweight: the operations and the Wyckoff positions of a setting are
// derived once and shared thereafter, so the first query of each setting pays
// for it. Priming moves that cost off the query path.
//
// A bitmask rather than a struct of flags: a caller names the families it
// wants, and the set is one value it can pass on. PointGroup and RodGroup are
// deliberately absent — they are built by value, not cached, so there is
// nothing to prime.
enum class Warm : unsigned {
  space_groups = 1U << 0, // the 530 Hall settings of the 230 space groups
  layer_groups = 1U << 1, // the 116 Hall settings of the 80 layer groups
  all = space_groups | layer_groups,
};

[[nodiscard]] constexpr Warm operator|(Warm a, Warm b) noexcept {
  return static_cast<Warm>(std::to_underlying(a) | std::to_underlying(b));
}

[[nodiscard]] constexpr bool contains(Warm set, Warm what) noexcept {
  return (std::to_underlying(set) & std::to_underlying(what)) != 0;
}

// Build the requested caches, one thread per family, and return once they are
// ready. Safe to call concurrently with any query: priming races the ordinary
// first-use path and both end up with the same shared objects.
void warmup(Warm what = Warm::all);

// The same, off the calling thread. The returned future must be waited on (or
// discarded deliberately) — priming continues either way.
[[nodiscard]] std::future<void> warmup_async(Warm what = Warm::all);

} // namespace cppcrystal
