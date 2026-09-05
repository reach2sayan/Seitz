#pragma once

#include "core/position_index.hpp"
#include "core/testable.hpp"
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

#include <boost/container/small_vector.hpp>

#include <cstdint>
#include <vector>

namespace cppcrystal {

// Tests candidate symmetry operations against a cell. Built once per cell and
// tolerance: the atoms are re-ordered by (type, distance to the nearest
// lattice point) so that the cheap rejection probes the most discriminating
// atoms first, and a PositionIndex over them makes the full check
// O(n log n) per operation instead of O(n^2).
//
// NOT thread-safe: check_total_overlap runs once per candidate (rotation,
// translation) pair, so its working buffers live here and are reused rather
// than reallocated per call. They are mutable because every caller holds the
// checker by const reference; the price is that one checker serves one thread.
class CPPCRYSTAL_TESTABLE OverlapChecker {
public:
  OverlapChecker(Cell const &cell, double symprec);

  // True iff x -> rot . x + trans maps the cell onto itself within symprec.
  [[nodiscard]] bool check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot) const;

private:
  // Cheap rejection: a few atoms must map onto some atom of their type. Takes
  // the transposed rotation rather than the mapped positions, so it can build
  // just the probe rows -- the full n x 3 image is only worth computing once
  // the probes have passed.
  [[nodiscard]] bool possible_overlap(Matrix3d const &rot_transposed,
                                      Vector3d const &trans) const;

  Cell sorted_; // the input cell with atoms sorted as described above
  PositionIndex index_; // over sorted_; ties the checker to its address

  mutable Positions rotated_; // n x 3 images
  mutable std::vector<boost::container::small_vector<int, 2>> images_;
  mutable std::vector<std::uint8_t> taken_;
  mutable PositionIndex::Scratch scratch_;
};

} // namespace cppcrystal
