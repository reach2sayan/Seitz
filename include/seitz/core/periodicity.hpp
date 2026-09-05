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

// 3D space group : {periodic,  periodic,  periodic}  (dim 3)
// layer group    : {periodic,  periodic,  aperiodic} (dim 2, aperiodic c)
// rod group      : {aperiodic, aperiodic, periodic}  (dim 1, periodic c)
// point group    : {aperiodic, aperiodic, aperiodic} (dim 0, a cluster)
using CellPeriodicity = std::array<AxisKind, 3>;

// The fully-periodic descriptor (a 3D space-group cell).
[[nodiscard]] constexpr CellPeriodicity all_periodic() noexcept {
  return {AxisKind::periodic, AxisKind::periodic, AxisKind::periodic};
}

// The family a cell's periodicity puts it in: any aperiodic axis = layer,
// all-periodic = 3D space group. (Rod and cluster cells belong to generation,
// not to the determination.)
[[nodiscard]] constexpr GroupFamily
family_of(CellPeriodicity const &p) noexcept {
  return std::ranges::contains(p, AxisKind::aperiodic) ? GroupFamily::layer
                                                       : GroupFamily::space;
}

// Layer group: periodic in the plane, aperiodic along `axis` (c in the
// conventional setting). Inverse of aperiodic_axis().
[[nodiscard]] constexpr CellPeriodicity aperiodic_along(int axis) noexcept {
  CellPeriodicity p = all_periodic();
  p[static_cast<std::size_t>(axis)] = AxisKind::aperiodic;
  return p;
}

// The cluster descriptor: no periodic axis at all (a 0D point set).
[[nodiscard]] constexpr CellPeriodicity none_periodic() noexcept {
  return {AxisKind::aperiodic, AxisKind::aperiodic, AxisKind::aperiodic};
}

// Rod group: periodic along `axis` only, the other two vacuum-padded. Dual of
// aperiodic_along().
[[nodiscard]] constexpr CellPeriodicity periodic_along(int axis) noexcept {
  CellPeriodicity p = none_periodic();
  p[static_cast<std::size_t>(axis)] = AxisKind::periodic;
  return p;
}

// The single aperiodic axis if there is exactly one — a layer group's c.
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

// The mixed-periodicity forms of the two folds below, out of line: their
// per-axis loop pushed the inline body past the inliner's budget, and these run
// once per atom per candidate operation, so the fully-periodic fast path has to
// fold into its caller.
[[nodiscard]] Vector3d minimal_image_mixed(Vector3d const &diff,
                                           CellPeriodicity const &p) noexcept;
[[nodiscard]] Vector3d wrap_mixed(Vector3d const &v,
                                  CellPeriodicity const &p) noexcept;

} // namespace detail

// Minimal-image fractional offset of `diff`: periodic components folded to
// their nearest-lattice-point residue, aperiodic ones left raw (no images).
[[nodiscard]] inline Vector3d minimal_image(Vector3d const &diff,
                                            CellPeriodicity const &p) noexcept {
  return p == all_periodic() ? Vector3d{math::nearest_offset(diff)}
                             : detail::minimal_image_mixed(diff, p);
}

// Fold into [0, 1) along every periodic axis, leaving aperiodic axes raw.
[[nodiscard]] inline Vector3d wrap(Vector3d const &v,
                                   CellPeriodicity const &p) noexcept {
  return p == all_periodic() ? math::wrap_to_unit_cell(v)
                             : detail::wrap_mixed(v, p);
}

} // namespace seitz

#pragma GCC visibility pop
