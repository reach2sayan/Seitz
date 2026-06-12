#pragma once

#include <array>
#include <optional>

namespace spglib {

// Whether a single lattice axis is periodic (translations along it are lattice
// translations) or aperiodic (it is a finite / vacuum direction). Replaces the
// magic "which axis is special" integer with a named two-state kind.
enum class AxisKind { periodic, aperiodic };

// The periodicity of all three axes at once — the general descriptor that covers
// every dimensionality the generation layer produces:
//   * 3D space group : {periodic,  periodic,  periodic}  (dim 3)
//   * layer group    : {periodic,  periodic,  aperiodic} (dim 2, aperiodic c)
//   * rod group      : {aperiodic, aperiodic, periodic}  (dim 1, periodic c)
//   * point group    : {aperiodic, aperiodic, aperiodic} (dim 0, a cluster)
// This supersedes Cell's single std::optional<int> aperiodic_axis, which can name
// at most one special axis and so cannot describe a rod (two aperiodic axes). The
// two bridge helpers below interoperate with the existing single-axis API while
// the wider migration is staged.
using CellPeriodicity = std::array<AxisKind, 3>;

// Number of periodic axes — the structure's periodic dimensionality (0..3).
[[nodiscard]] constexpr int
periodic_dimension(CellPeriodicity const &p) noexcept {
  int dim = 0;
  for (AxisKind const k : p) {
    if (k == AxisKind::periodic) {
      ++dim;
    }
  }
  return dim;
}

// The fully-periodic descriptor (a 3D space-group cell).
[[nodiscard]] constexpr CellPeriodicity all_periodic() noexcept {
  return {AxisKind::periodic, AxisKind::periodic, AxisKind::periodic};
}

// Build a periodicity from the legacy single-aperiodic-axis convention: nullopt
// (3D, all periodic) or one aperiodic axis (a layer group's c). Lets the new
// type interoperate with code still passing std::optional<int>.
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
  for (int axis = 0; axis < 3; ++axis) {
    if (p[static_cast<std::size_t>(axis)] == AxisKind::aperiodic) {
      if (found) {
        return std::nullopt; // more than one aperiodic axis (rod / point)
      }
      found = axis;
    }
  }
  return found;
}

} // namespace spglib
