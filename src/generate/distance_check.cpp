#include <cppcrystal/generate/distance_check.hpp>

#include "core/position_index.hpp"
#include "math/integer_matrix.hpp"
#include <cppcrystal/data/element_data.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>

namespace cppcrystal::generate {

namespace {

// Covalent radius of an atom type, falling back to tol.fallback_radius for an
// untabulated element so the check stays well-defined. constexpr — a pure
// lookup over the compile-time covalent-radius table.
[[nodiscard]] constexpr double
radius_of(int type, DistanceTolerance const &tol) noexcept {
  return data::covalent_radius(type).value_or(tol.fallback_radius);
}

// The shared check of the validity tests: every pair of distinct atoms must be
// at least scale * (r_i + r_j) apart under the family's metric, and every atom
// at least 2 * scale * r_i from its own periodic images. No threshold exceeds
// the cutoff 2 * scale * max radius, so a bucket grid with that edge yields
// every pair that could possibly be too close and nothing else is measured.
[[nodiscard]] bool pairwise_distances_ok(Positions const &positions,
                                         Types const &types,
                                         Matrix3d const &lattice,
                                         CellPeriodicity const &periodicity,
                                         DistanceTolerance tol) noexcept {
  if (types.empty()) {
    return true;
  }
  auto const radius = types | std::views::transform([&](int type) {
                        return radius_of(type, tol);
                      });
  double const cutoff = 2.0 * tol.scale * std::ranges::max(radius);

  // Self-images: the nearest non-trivial image of any point is the same
  // lattice vector for every atom (infinite with no periodic axis).
  double const self_image =
      minimum_image_distance(Vector3d::Zero(), Vector3d::Zero(), lattice,
                             periodicity, Images::nontrivial);
  if (std::ranges::any_of(
          radius, [&](double r) { return self_image < 2.0 * tol.scale * r; })) {
    return false;
  }

  PositionIndex const index(positions, types, lattice, cutoff, periodicity);
  auto const row = [&](Index i) { return Vector3d(positions.row(i)); };
  for (Index i = 0; i < positions.rows(); ++i) {
    auto later = index.candidates(row(i)) |
                 std::views::filter([i](int j) { return j > i; });
    for (int const j : later) {
      double const min_dist = tol.scale * (radius[i] + radius[j]);
      if (minimum_image_distance(row(i), row(j), lattice, periodicity) <
          min_dist) {
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
                              Images images) noexcept {
  // Fold each periodic component into [-0.5, 0.5]; keep aperiodic components
  // raw (no images along them).
  Vector3d const base = minimal_image(a - b, periodicity);

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
    if (images == Images::nontrivial && off.squaredNorm() < 1e-20) {
      continue; // the home image of an atom against itself
    }
    best = std::min(best, (lattice * off).squaredNorm());
  }
  return std::sqrt(best);
}

bool distances_valid(Cell const &cell, DistanceTolerance tol) noexcept {
  return pairwise_distances_ok(cell.positions(), cell.types(),
                               cell.lattice().matrix(), cell.periodicity(),
                               tol);
}

} // namespace cppcrystal::generate
