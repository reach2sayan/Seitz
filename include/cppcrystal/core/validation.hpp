#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/tolerance.hpp>

#include <cmath>

namespace cppcrystal {

// Reject degenerate cells before the symmetry pipeline touches them: an empty
// cell, or a (near-)singular lattice whose inverse would otherwise propagate
// NaN/inf into the public result. One cheap up-front guard at each entry point
// keeps the internal Matrix3d::inverse() calls safe without per-call checks.
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

} // namespace cppcrystal
