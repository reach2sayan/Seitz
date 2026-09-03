#include <cppcrystal/group/wyckoff_position.hpp>

#include "core/matrix_order.hpp"
#include <cppcrystal/core/tolerance.hpp>
#include "math/fractional.hpp"

#include <vector>

namespace cppcrystal::group {

char WyckoffPosition::letter() const noexcept {
  return letter_ < 26 ? static_cast<char>('a' + letter_)
                      : static_cast<char>('A' + (letter_ - 26));
}

Positions WyckoffPosition::get_all_positions(Vector3d const &xyz) const {
  // Project the input onto the position's canonical coordinate so a point only
  // approximately on the position still yields the exact orbit; exact
  // (rational) database coordinates are compared at the default symprec.
  Vector3d const canonical = repr_rotation_.cast<double>() * xyz +
                             repr_translation_;

  std::vector<Vector3d> orbit;
  orbit.reserve(orbit_ops_.size());
  for (auto const &op : orbit_ops_) {
    push_unique(orbit, math::wrap_to_unit_cell(op.apply(canonical)),
                [](Vector3d const &p, Vector3d const &q) {
                  return math::same_point(p, q, kDefaultSymprec);
                });
  }

  return to_positions(orbit);
}

} // namespace cppcrystal::group
