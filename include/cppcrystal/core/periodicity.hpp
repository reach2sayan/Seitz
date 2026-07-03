#pragma once

#include <algorithm>
#include <array>
#include <optional>
#include <ranges>

namespace cppcrystal {

enum class AxisKind { periodic, aperiodic };

//   3D space group : {periodic,  periodic,  periodic}  (dim 3)
//   layer group    : {periodic,  periodic,  aperiodic} (dim 2, aperiodic c)
//   rod group      : {aperiodic, aperiodic, periodic}  (dim 1, periodic c)
//   point group    : {aperiodic, aperiodic, aperiodic} (dim 0, a cluster)
using CellPeriodicity = std::array<AxisKind, 3>;

[[nodiscard]] constexpr int
periodic_dimension(CellPeriodicity const &p) noexcept {
  return static_cast<int>(std::ranges::count(p, AxisKind::periodic));
}

// The fully-periodic descriptor (a 3D space-group cell).
[[nodiscard]] constexpr CellPeriodicity all_periodic() noexcept {
  return {AxisKind::periodic, AxisKind::periodic, AxisKind::periodic};
}

[[nodiscard]] constexpr CellPeriodicity
periodicity_from_aperiodic_axis(std::optional<int> aperiodic_axis) noexcept {
  CellPeriodicity p = all_periodic();
  if (aperiodic_axis) {
    p[static_cast<std::size_t>(*aperiodic_axis)] = AxisKind::aperiodic;
  }
  return p;
}

// The inverse bridge: the single aperiodic axis if there is exactly one,
// otherwise std::nullopt (all-periodic 3D, or the >1-aperiodic rod/point cases
// the legacy API cannot represent). A caller that must stay on the single-axis
// API can detect those cases by the nullopt-with-non-3D mismatch.
[[nodiscard]] constexpr std::optional<int>
single_aperiodic_axis(CellPeriodicity const &p) noexcept {
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

} // namespace cppcrystal
