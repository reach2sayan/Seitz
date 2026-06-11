#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/types.hpp>

#include <optional>
#include <vector>

namespace spglib {

// True when fractional positions a and b coincide modulo the lattice, measured
// by Cartesian minimal-image distance <= symprec (spglib cel_is_overlap /
// has_overlap). With `aperiodic_axis` set (layer groups), that axis is NOT
// periodic: only the two periodic components are folded to the minimal image;
// the aperiodic component keeps its raw difference (cel_layer_is_overlap).
[[nodiscard]] bool is_overlap(Vector3d const &a, Vector3d const &b,
                              Matrix3d const &lattice, double symprec,
                              std::optional<int> aperiodic_axis =
                                  std::nullopt) noexcept;

[[nodiscard]] inline bool
is_overlap_same_type(Vector3d const &a, Vector3d const &b, int type_a,
                     int type_b, Matrix3d const &lattice, double symprec,
                     std::optional<int> aperiodic_axis = std::nullopt) noexcept {
  return type_a == type_b &&
         is_overlap(a, b, lattice, symprec, aperiodic_axis);
}

// Precomputes a distance-sorted copy of a cell so that many candidate symmetry
// operations can be tested efficiently against it. Port of overlap.c's
// OverlapChecker, minus the layer-group path (space groups only).
class OverlapChecker {
public:
  explicit OverlapChecker(Cell const &cell);

  // True iff x -> rot . x + trans maps the cell onto itself within symprec.
  // (The identity-rotation fast path is detected internally.)
  [[nodiscard]] bool check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot,
                                         double symprec) const;

private:
  // Cheap rejection: confirm a few atoms map onto some atom before the full
  // check.
  [[nodiscard]] bool possible_overlap(Vector3d const &trans,
                                      Matrix3i const &rot,
                                      double symprec) const;

  Matrix3d lattice_;
  Positions pos_sorted_; // sorted by distance to nearest lattice point
  std::vector<int> types_sorted_;
  int size_;
  std::optional<int> aperiodic_axis_; // set for layer cells; nullopt for 3D
};

} // namespace spglib
