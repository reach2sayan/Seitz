#pragma once

#include <seitz/core/error.hpp>
#include <seitz/group/wyckoff.hpp>

#pragma GCC visibility push(default)

namespace seitz::group {

// A rod group (1..75) as a standalone, structure-free object: operations
// decoded from the generated rod tables, Wyckoff positions derived from them.
// The 1D-periodic sibling of SpaceGroup (3D), SpaceGroup::from_layer_* (2D)
// and PointGroup (0D).
//
// Positions (general and special) come from the affine fixed-locus
// arrangement: each op's locus solves (R - I) p = -t + n e_axis modulo the rod
// lattice, the loci are closed under intersection, one position per orbit of
// loci.
class RodGroup : public GroupBase {
public:
  // By rod-group number (1..75), over the rod tables transcribed from
  // tools/rod_group_data.csv.
  [[nodiscard]] static Result<RodGroup> from_number(int number);

  // The one periodic axis (0=a, 1=b, 2=c); orbit expansion folds only this
  // component, the other two being vacuum directions.
  [[nodiscard]] int periodic_axis() const noexcept { return periodic_axis_; }

private:
  RodGroup() = default;
  int periodic_axis_ = 2;
};

} // namespace seitz::group

#pragma GCC visibility pop
