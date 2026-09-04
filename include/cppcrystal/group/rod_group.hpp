#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/group/wyckoff.hpp>

// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal::group {

// A rod group (1..75) as a standalone, structure-free object: it owns its
// symmetry operations (decoded from the generated rod tables) and its Wyckoff
// positions, derived from those operations. The 1D-periodic sibling of
// group::SpaceGroup (3D), SpaceGroup::from_layer_* (2D) and group::PointGroup
// (0D).
//
// The Wyckoff positions (general AND special) are derived from the operations
// by the affine fixed-locus arrangement: with the periodic axis handled modulo
// the rod lattice, each operation's fixed locus is the affine solution of
// (R - I) p = -t + n*e_axis, the loci are closed under intersection, and one
// Wyckoff position is emitted per orbit of loci.
class RodGroup : public GroupBase {
public:
  // Build by rod-group number (1..75). Requires the generated rod tables
  // (src/data/rod_group_tables.hpp via tools/transcribe_rod_groups.py).
  [[nodiscard]] static Result<RodGroup> from_number(int number);

  // The single periodic axis (0=a, 1=b, 2=c); orbit expansion folds only this
  // component (the other two are aperiodic / vacuum directions).
  [[nodiscard]] int periodic_axis() const noexcept { return periodic_axis_; }

private:
  RodGroup() = default;
  int periodic_axis_ = 2;
};

} // namespace cppcrystal::group

#pragma GCC visibility pop
