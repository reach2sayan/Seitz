#pragma once

#include <seitz/core/error.hpp>
#include <seitz/group/wyckoff.hpp>

#include <string_view>

#pragma GCC visibility push(default)

namespace seitz::group {

// A crystallographic point group (1..32) as a standalone, structure-free
// object:
// Owns symmetry operations (pure rotations / rotoinversions about the origin,
// expressed as integer matrices in the conventional basis)
// + Wyckoff positions.
//
// Operations : a representative symmorphic space group of the point group;
// then Wyckoff positions : arrangement of fixed subspaces
// (the origin, the rotation axes, the mirror planes, and their intersections).
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

} // namespace seitz::group

#pragma GCC visibility pop
