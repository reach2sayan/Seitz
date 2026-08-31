#include <cppcrystal/generate/distance_check.hpp>

#include <cppcrystal/data/element_data.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <algorithm>
#include <limits>
#include <ranges>
#include <vector>

namespace cppcrystal::generate {

namespace {

// Covalent radius of an atom type, falling back to tol.fallback_radius for an
// untabulated element so the check stays well-defined. constexpr — a pure
// lookup over the compile-time covalent-radius table.
[[nodiscard]] constexpr double radius_of(int type,
                                         DistanceTolerance const &tol) noexcept {
  return data::covalent_radius(type).value_or(tol.fallback_radius);
}

// The shared pair loop of the validity checks: j == i covers an atom against
// its own periodic images (cell-too-small); j > i covers distinct pairs and
// their images. `distance(i, j)` is the family's metric.
template <class Metric>
[[nodiscard]] bool pairwise_distances_ok(Types const &types,
                                         DistanceTolerance tol,
                                         Metric &&distance) noexcept {
  Index const n = static_cast<Index>(types.size());
  std::vector<double> radius(types.size());
  std::ranges::transform(types, radius.begin(),
                         [&](int type) { return radius_of(type, tol); });

  for (Index i = 0; i < n; ++i) {
    double const ri = radius[static_cast<std::size_t>(i)];
    for (Index j = i; j < n; ++j) {
      double const min_dist =
          tol.scale * (ri + radius[static_cast<std::size_t>(j)]);
      if (distance(i, j) < min_dist) {
        return false;
      }
    }
  }
  return true;
}

} // namespace

double minimum_image_distance(Vector3d const &a, Vector3d const &b,
                              Matrix3d const &lattice,
                              CellPeriodicity const &periodicity,
                              bool include_origin) noexcept {
  // Fold each periodic component into [-0.5, 0.5]; keep aperiodic components
  // raw (no images along them).
  Vector3d const diff = a - b;
  Vector3d base = diff;
  for (auto const [axis, kind] : periodicity | std::views::enumerate) {
    if (kind == AxisKind::periodic) {
      base[axis] = diff[axis] - math::nint(diff[axis]);
    }
  }

  // Search the adjacent cells along the periodic axes only: for skewed
  // lattices the true nearest image can sit one cell over.
  auto const neighbours = [&](std::size_t axis) {
    return periodicity[axis] == AxisKind::periodic ? std::views::iota(-1, 2)
                                                   : std::views::iota(0, 1);
  };
  double best = std::numeric_limits<double>::infinity();
  for (auto const [n0, n1, n2] : std::views::cartesian_product(
           neighbours(0), neighbours(1), neighbours(2))) {
    Vector3d const off =
        base + Vector3d(static_cast<double>(n0), static_cast<double>(n1),
                        static_cast<double>(n2));
    if (!include_origin && off.squaredNorm() < 1e-20) {
      continue; // the home image of an atom against itself
    }
    best = std::min(best, (lattice * off).squaredNorm());
  }
  return std::sqrt(best);
}

bool distances_valid(Cell const &cell, DistanceTolerance tol) noexcept {
  return pairwise_distances_ok(cell.types(), tol, [&](Index i, Index j) {
    return minimum_image_distance(cell.position(i), cell.position(j),
                                  cell.lattice(), cell.periodicity(),
                                  /*include_origin=*/i != j);
  });
}

bool cluster_distances_valid(Positions const &coordinates, Types const &types,
                             DistanceTolerance tol) noexcept {
  constexpr CellPeriodicity kNonPeriodic = {
      AxisKind::aperiodic, AxisKind::aperiodic, AxisKind::aperiodic};
  return pairwise_distances_ok(types, tol, [&](Index i, Index j) {
    return minimum_image_distance(coordinates.row(i).transpose(),
                                  coordinates.row(j).transpose(),
                                  Matrix3d::Identity(), kNonPeriodic,
                                  /*include_origin=*/i != j);
  });
}

} // namespace cppcrystal::generate
