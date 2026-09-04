#include <cppcrystal/generate/distance_check.hpp>

#include "core/position_index.hpp"
#include "math/integer_matrix.hpp"
#include <cppcrystal/data/element_data.hpp>

#include <boost/container/static_vector.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

namespace cppcrystal::generate {

namespace {

// The neighbour images to search, built once per cell: each entry pairs a
// fractional offset n with its Cartesian image lattice * n. Because
// lattice * (base + n) == lattice * base + lattice * n, the per-pair loop then
// costs one vector add and a squaredNorm per image instead of a fresh 3x3
// matvec -- and the pair loop runs once per atom pair per attempt.
//
// The distributed form is not bit-identical to the matvec it replaces (the two
// round differently in the last place). That is acceptable here and only here:
// the comparison this feeds is `distance < scale * (r_i + r_j)`, a threshold
// from a covalent-radius table with two significant figures, not a
// tolerance-parity check against the reference implementation. The seeded
// generation tests pin the consequence -- a structure accepted before must
// still be accepted, for the same seed.
class MinimumImage {
public:
  MinimumImage(Matrix3d const &lattice, CellPeriodicity const &periodicity)
      : lattice_{lattice}, periodicity_{periodicity} {
    // Adjacent cells along the periodic axes only: for skewed lattices the
    // true nearest image can sit one cell over, while an aperiodic axis has no
    // images at all.
    auto const neighbours = [&](std::size_t axis) {
      return periodicity[axis] == AxisKind::periodic ? std::views::iota(-1, 2)
                                                     : std::views::iota(0, 1);
    };
    for (auto const [n0, n1, n2] : std::views::cartesian_product(
             neighbours(0), neighbours(1), neighbours(2))) {
      Vector3d const offset(static_cast<double>(n0), static_cast<double>(n1),
                            static_cast<double>(n2));
      shifts_.emplace_back(offset, Vector3d(lattice * offset));
    }
  }

  [[nodiscard]] double distance(Vector3d const &a, Vector3d const &b,
                                Images images) const noexcept {
    // Fold each periodic component into [-0.5, 0.5]; keep aperiodic components
    // raw (no images along them).
    Vector3d const base = minimal_image(a - b, periodicity_);
    Vector3d const home = lattice_ * base;
    double best = std::numeric_limits<double>::infinity();
    for (auto const &[offset, shift] : shifts_) {
      if (images == Images::nontrivial &&
          (base + offset).squaredNorm() < 1e-20) {
        continue; // the home image of an atom against itself
      }
      best = std::min(best, (home + shift).squaredNorm());
    }
    return std::sqrt(best);
  }

private:
  Matrix3d lattice_;
  CellPeriodicity periodicity_;
  boost::container::static_vector<std::pair<Vector3d, Vector3d>, 27> shifts_;
};

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
                                         DistanceTolerance tol) {
  if (types.empty()) {
    return true;
  }
  // Materialised, not a lazy transform view: radius[i] / radius[j] are read in
  // the pair loop below, and a view would re-run the covalent-radius lookup on
  // every one of those reads.
  std::vector<double> const radius{std::from_range,
                                   types | std::views::transform([&](int type) {
                                     return radius_of(type, tol);
                                   })};
  double const cutoff = 2.0 * tol.scale * std::ranges::max(radius);

  MinimumImage const metric(lattice, periodicity);

  // Self-images: the nearest non-trivial image of any point is the same
  // lattice vector for every atom (infinite with no periodic axis).
  double const self_image = metric.distance(Vector3d::Zero(), Vector3d::Zero(),
                                            Images::nontrivial);
  if (std::ranges::any_of(
          radius, [&](double r) { return self_image < 2.0 * tol.scale * r; })) {
    return false;
  }

  PositionIndex const index(positions, types, lattice, cutoff, periodicity);
  PositionIndex::Scratch scratch;
  auto const row = [&](Index i) { return Vector3d(positions.row(i)); };
  for (Index i = 0; i < positions.rows(); ++i) {
    Vector3d const from = row(i);
    auto later = index.candidates(from, scratch) |
                 std::views::filter([i](int j) { return j > i; });
    for (int const j : later) {
      double const min_dist =
          tol.scale * (radius[static_cast<std::size_t>(i)] +
                       radius[static_cast<std::size_t>(j)]);
      if (metric.distance(from, row(j), Images::all) < min_dist) {
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
  return MinimumImage{lattice, periodicity}.distance(a, b, images);
}

bool distances_valid(Cell const &cell, DistanceTolerance tol) {
  return pairwise_distances_ok(cell.positions(), cell.types(),
                               cell.lattice().matrix(), cell.periodicity(),
                               tol);
}

} // namespace cppcrystal::generate
