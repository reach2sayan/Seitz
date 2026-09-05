#pragma once

#include <seitz/alloy/cluster.hpp>
#include <seitz/core/fractional.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <algorithm>
#include <ranges>

// Private to src/alloy: the three geometric primitives the cluster, pool and
// CVM builders share. All fractional, all on a caller-supplied tolerance --
// each is a line over seitz::math or a view.
namespace seitz::alloy::detail {

// A fractional delta that is a lattice translation: every component an integer
// to within `tol`. In the parent's own fractional frame this is what "modulo
// the cell" means, so it is a rounding test and not a metric one.
[[nodiscard]] inline bool integral(Vector3d const &delta, double tol) noexcept {
  return approx_zero(math::nearest_offset(delta), tol);
}

// The point of `c` at `position` -- exactly, with no modulo, because every
// caller has already applied the whole-cluster shift -- or `c.end()`. An
// iterator rather than an index: the callers want the point's own fields, and
// the two that do want a position recover it with ranges::distance.
[[nodiscard]] inline auto coincident(Cluster const &c, Vector3d const &position,
                                     double tol) {
  return std::ranges::find_if(c, [&position, tol](ClusterPoint const &p) {
    return (p.position - position).norm() < tol;
  });
}

// The integer translations of the inclusive box [lo, hi], LAST component
// varying fastest. The walk order is load-bearing: the correlation sums
// downstream accumulate in it, so reordering this changes results in the last
// bits even though the set of translations is the same.
[[nodiscard]] inline auto lattice_box(Vector3i const &lo,
                                      Vector3i const &hi) noexcept {
  return std::views::cartesian_product(std::views::iota(lo(0), hi(0) + 1),
                                       std::views::iota(lo(1), hi(1) + 1),
                                       std::views::iota(lo(2), hi(2) + 1));
}

[[nodiscard]] inline auto lattice_box(int reach) noexcept {
  return lattice_box(Vector3i::Constant(-reach), Vector3i::Constant(reach));
}

// The mixed-radix bound of a cluster's points under `bound`: the shape of the
// odometer that walks its decorations, its configurations or its subsets.
// Materialized because detail::mixed_radix takes a contiguous span.
template <class Bound>
  requires std::convertible_to<
      std::invoke_result_t<Bound &,
                           std::ranges::range_reference_t<Cluster const &>>,
      int>
[[nodiscard]] std::vector<int> radix_of(Cluster const &c, Bound bound) {
  return std::vector<int>{std::from_range,
                          c | std::views::transform(bound)};
}

} // namespace seitz::alloy::detail
