#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>

#include <cmath>

namespace spglib {

// Reject degenerate cells before the symmetry pipeline touches them (Phase 9
// hardening): an empty cell, or a (near-)singular lattice whose inverse would
// otherwise propagate NaN/inf into the public result. One cheap up-front guard
// at each entry point keeps the determination's internal Matrix3d::inverse()
// calls safe without per-call checks.
[[nodiscard]] inline Result<void> validate_cell(Cell const &cell) {
  if (cell.size() == 0) {
    return leaf::new_error(e_empty_cell{});
  }
  double const determinant = cell.lattice().determinant();
  if (std::abs(determinant) < kZeroPrec) {
    return leaf::new_error(e_invalid_lattice{determinant});
  }
  return {};
}

} // namespace spglib
