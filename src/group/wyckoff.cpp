#include <cppcrystal/group/wyckoff.hpp>

#include <cppcrystal/core/tolerance.hpp>

#include <algorithm>
#include <ranges>

namespace cppcrystal::group {

Positions Wyckoff::orbit(Vector3d const &xyz,
                         CellPeriodicity const &periodicity) const {
  Vector3d const seed = canonical(xyz);

  // Written straight into the result rather than into a vector<Vector3d> that
  // to_positions would then copy row by row: there is at most one point per
  // coset representative, and for a generic `xyz` -- the case the generator
  // asks for, once per placement per attempt -- it is exactly one each, so the
  // shrink below never runs.
  Positions points(static_cast<Index>(orbit_ops_.size()), 3);
  Index kept = 0;
  for (auto const &op : orbit_ops_) {
    Vector3d const point = wrap(op.apply(seed), periodicity);
    auto const coincides = [&](Index k) {
      return minimal_image(Vector3d(points.row(k).transpose() - point),
                           periodicity)
                 .cwiseAbs()
                 .maxCoeff() < kDefaultSymprec;
    };
    if (!std::ranges::any_of(std::views::iota(Index{0}, kept), coincides)) {
      points.row(kept++) = point.transpose();
    }
  }
  if (kept < points.rows()) {
    points.conservativeResize(kept, 3); // the accidental-coincidence case
  }
  return points;
}

} // namespace cppcrystal::group
