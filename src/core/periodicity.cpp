#include <cppcrystal/core/periodicity.hpp>

#include <ranges>

namespace cppcrystal::detail {

Vector3d minimal_image_mixed(Vector3d const &diff,
                             CellPeriodicity const &p) noexcept {
  Vector3d out = math::nearest_offset(diff);
  for (auto const [axis, kind] : p | std::views::enumerate) {
    if (kind == AxisKind::aperiodic) {
      out[axis] = diff[axis];
    }
  }
  return out;
}

Vector3d wrap_mixed(Vector3d const &v, CellPeriodicity const &p) noexcept {
  Vector3d out;
  for (auto const [axis, kind] : p | std::views::enumerate) {
    out[axis] = kind == AxisKind::aperiodic ? v[axis]
                                            : math::wrap_to_unit_cell(v[axis]);
  }
  return out;
}

} // namespace cppcrystal::detail
