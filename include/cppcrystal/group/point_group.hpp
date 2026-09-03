#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/group/wyckoff.hpp>

#include <string_view>

namespace cppcrystal::group {

// A crystallographic point group (1..32) as a standalone, structure-free
// object: it owns its symmetry operations (pure rotations / rotoinversions
// about the origin, expressed as integer matrices in the conventional basis)
// and its Wyckoff positions as first-class objects.
//
// The operations are taken from a representative symmorphic space group of the
// point group (the distinct rotation parts of its conventional operations); the
// Wyckoff positions are then derived from the arrangement of fixed subspaces
// (the origin, the rotation axes, the mirror planes, and their intersections).
// This is the 0D counterpart of SpaceGroup (3D) and SpaceGroup::from_layer_*
// (2D) for crystal/cluster generation.
class PointGroup : public GroupBase {
public:
  // Build by international point-group number (1..32, the spglib pointgroup
  // numbering: 1, -1, 2, m, 2/m, 222, ...).
  [[nodiscard]] static Result<PointGroup> from_number(int number);

  [[nodiscard]] std::string_view schoenflies() const noexcept {
    return schoenflies_;
  }

  // A space group whose point group is this one (a symmorphic P group); its
  // crystal system fixes the metric used when realising a cluster in Cartesian
  // space (generate::random_cluster).
  [[nodiscard]] int representative_spacegroup() const noexcept {
    return representative_spacegroup_;
  }

private:
  PointGroup() = default;

  std::string_view schoenflies_;
  int representative_spacegroup_ = 0;
};

} // namespace cppcrystal::group
