#pragma once

#include "core/position_index.hpp"
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

namespace cppcrystal {

// Tests candidate symmetry operations against a cell. Built once per cell and
// tolerance: the atoms are re-ordered by (type, distance to the nearest
// lattice point) so that the cheap rejection probes the most discriminating
// atoms first, and a PositionIndex over them makes the full check
// O(n log n) per operation instead of O(n^2).
class OverlapChecker {
public:
  OverlapChecker(Cell const &cell, double symprec);

  // True iff x -> rot . x + trans maps the cell onto itself within symprec.
  [[nodiscard]] bool check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot) const;

private:
  // Cheap rejection: a few atoms must map onto some atom of their type.
  [[nodiscard]] bool possible_overlap(Positions const &rotated) const;

  Cell sorted_; // the input cell with atoms sorted as described above
  double symprec_;
  PositionIndex index_; // over sorted_; ties the checker to its address
};

} // namespace cppcrystal
