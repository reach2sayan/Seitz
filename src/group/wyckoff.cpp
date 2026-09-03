#include <cppcrystal/group/wyckoff.hpp>

#include <cppcrystal/core/tolerance.hpp>

#include "core/matrix_order.hpp"

#include <vector>

namespace cppcrystal::group {

Positions Wyckoff::orbit(Vector3d const &xyz,
                         CellPeriodicity const &periodicity) const {
  Vector3d const seed = canonical(xyz);
  std::vector<Vector3d> points;
  points.reserve(orbit_ops_.size());
  for (auto const &op : orbit_ops_) {
    push_unique(points, wrap(op.apply(seed), periodicity),
                [&](Vector3d const &p, Vector3d const &q) {
                  return minimal_image(Vector3d(p - q), periodicity)
                             .cwiseAbs()
                             .maxCoeff() < kDefaultSymprec;
                });
  }
  return to_positions(points);
}

} // namespace cppcrystal::group
