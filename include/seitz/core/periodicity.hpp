#pragma once

#include <seitz/core/fractional.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/types.hpp>

#include <algorithm>
#include <array>
#include <optional>
#include <ranges>

#pragma GCC visibility push(default)

namespace seitz {

enum class AxisKind { periodic, aperiodic };

//   3D space group : {periodic,  periodic,  periodic}  (dim 3)
//   layer group    : {periodic,  periodic,  aperiodic} (dim 2, aperiodic c)
//   rod group      : {aperiodic, aperiodic, periodic}  (dim 1, periodic c)
//   point group    : {aperiodic, aperiodic, aperiodic} (dim 0, a cluster)
using CellPeriodicity = std::array<AxisKind, 3>;

// The fully-periodic descriptor (a 3D space-group cell).
[[nodiscard]] constexpr CellPeriodicity all_periodic() noexcept {
  return {AxisKind::periodic, AxisKind::periodic, AxisKind::periodic};
}

// The family a cell's periodicity puts it in: one aperiodic axis is a layer
// group, all-periodic is a 3D space group. (Rod and cluster cells are handled
// by the generation layer, not the space-group search.)
[[nodiscard]] constexpr GroupFamily
family_of(CellPeriodicity const &p) noexcept {
  return std::ranges::contains(p, AxisKind::aperiodic) ? GroupFamily::layer
                                                       : GroupFamily::space;
}

// The layer-group descriptor: periodic in the plane, aperiodic along `axis`
// (c, axis 2, in the conventional setting). The inverse of aperiodic_axis().
[[nodiscard]] constexpr CellPeriodicity aperiodic_along(int axis) noexcept {
  CellPeriodicity p = all_periodic();
  p[static_cast<std::size_t>(axis)] = AxisKind::aperiodic;
  return p;
}

// The cluster descriptor: no periodic axis at all (a 0D point set).
[[nodiscard]] constexpr CellPeriodicity none_periodic() noexcept {
  return {AxisKind::aperiodic, AxisKind::aperiodic, AxisKind::aperiodic};
}

// The rod descriptor: periodic along `axis` only, the two remaining axes being
// the vacuum-padded cross-section. The dual of aperiodic_along().
[[nodiscard]] constexpr CellPeriodicity periodic_along(int axis) noexcept {
  CellPeriodicity p = none_periodic();
  p[static_cast<std::size_t>(axis)] = AxisKind::periodic;
  return p;
}

// The single aperiodic axis if there is exactly one — a layer group's c. The
// layer path needs the axis index itself (to pick the in-plane pair, to reject
// the cubic point groups, ...); this is the one place that search lives.
// std::nullopt for the 3D case and for the rod/cluster cases, which have no
// single distinguished aperiodic axis.
[[nodiscard]] constexpr std::optional<int>
aperiodic_axis(CellPeriodicity const &p) noexcept {
  std::optional<int> found;

  for (auto const [axis, kind] : p | std::views::enumerate) {
    if (kind != AxisKind::aperiodic) {
      continue;
    }
    if (found) {
      return std::nullopt;
    }
    found = static_cast<int>(axis);
  }

  return found;
}

namespace detail {

// The mixed-periodicity forms of the two folds below, out of line. A layer,
// rod or cluster cell is the rare case, and its per-axis loop is what pushed
// the whole function past the inliner's budget: keeping it out of the header
// leaves an inline body small enough that the fully-periodic fast path folds
// into its caller, which matters because these run once per atom per candidate
// operation.
[[nodiscard]] Vector3d minimal_image_mixed(Vector3d const &diff,
                                           CellPeriodicity const &p) noexcept;
[[nodiscard]] Vector3d wrap_mixed(Vector3d const &v,
                                  CellPeriodicity const &p) noexcept;

} // namespace detail

// Minimal-image fractional offset of `diff`: every periodic component is
// folded to its nearest-lattice-point residue, every aperiodic component is
// left as the raw difference (no images along it).
[[nodiscard]] inline Vector3d minimal_image(Vector3d const &diff,
                                            CellPeriodicity const &p) noexcept {
  return p == all_periodic() ? Vector3d{math::nearest_offset(diff)}
                             : detail::minimal_image_mixed(diff, p);
}

// Fold a fractional coordinate into the cell [0, 1) along every periodic axis,
// leaving the aperiodic axes at their raw values — a layer/rod/cluster cell is
// not periodic along those.
[[nodiscard]] inline Vector3d wrap(Vector3d const &v,
                                   CellPeriodicity const &p) noexcept {
  return p == all_periodic() ? math::wrap_to_unit_cell(v)
                             : detail::wrap_mixed(v, p);
}

} // namespace seitz

#pragma GCC visibility pop
