#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/position_index.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>

namespace cppcrystal {

// True when fractional positions a and b coincide modulo the lattice, measured
// by Cartesian minimal-image distance <= symprec. With `aperiodic_axis` set
// (layer groups), that axis is NOT periodic: only the two periodic components
// are folded to the minimal image; the aperiodic component keeps its raw
// difference. The single-aperiodic-axis spelling of `coincident`.
[[nodiscard]] inline bool
is_overlap(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
           double symprec,
           std::optional<int> aperiodic_axis = std::nullopt) noexcept {
  return coincident(a, b, lattice, symprec,
                    periodicity_from_aperiodic_axis(aperiodic_axis));
}

[[nodiscard]] inline bool is_overlap_same_type(
    Vector3d const &a, Vector3d const &b, int type_a, int type_b,
    Matrix3d const &lattice, double symprec,
    std::optional<int> aperiodic_axis = std::nullopt) noexcept {
  return type_a == type_b && is_overlap(a, b, lattice, symprec, aperiodic_axis);
}

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
