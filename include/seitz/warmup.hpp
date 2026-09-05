#pragma once

#include <future>
#include <utility>

#pragma GCC visibility push(default)

namespace seitz {

// Which per-setting caches to build ahead of first use. group::SpaceGroup::of
// is a flyweight, so the first query of a setting pays for its operations and
// Wyckoff positions; priming moves that off the query path. A bitmask, so a
// caller passes the set of families on as one value. PointGroup and RodGroup
// are built by value, not cached, so there is nothing to prime.
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

// Build the requested caches, one thread per family, returning when ready.
// Safe to call concurrently with any query: priming and first use race to the
// same shared objects.
void warmup(Warm what = Warm::all);

// The same off the calling thread; the future must be waited on or discarded
// deliberately, priming continues either way.
[[nodiscard]] std::future<void> warmup_async(Warm what = Warm::all);

} // namespace seitz

#pragma GCC visibility pop
