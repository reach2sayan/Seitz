#include <seitz/generate/distance_check.hpp>

#include "core/position_index.hpp"
#include "math/integer_matrix.hpp"

#include <boost/container/flat_set.hpp>
#include <boost/container/static_vector.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

namespace seitz::generate {

namespace {

// The neighbour images to search, built once per cell: offset n paired with
// lattice * n. Since lattice * (base + n) = lattice * base + lattice * n, each
// image costs one add and a squaredNorm instead of a 3x3 matvec, in a loop that
// runs once per atom pair per attempt.
//
// The distributed form differs from the matvec in the last place. Acceptable
// here only: it feeds `distance < scale * (r_i + r_j)`, a two-significant-figure
// covalent-radius threshold, not a tolerance-parity check against the reference.
// The seeded generation tests pin the consequence.
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

// The shared check of the validity tests: every pair of distinct atoms must
// clear its type-pair minimum distance under the family's metric, and every
// atom the distance of its type against itself from its own periodic images.
// No threshold exceeds the largest minimum distance among the types present,
// so a bucket grid with that edge yields every pair that could possibly be too
// close and nothing else is measured.
[[nodiscard]] bool pairwise_distances_ok(Positions const &positions,
                                         Types const &types,
                                         Matrix3d const &lattice,
                                         CellPeriodicity const &periodicity,
                                         DistanceTolerance const &tol) {
  if (types.empty()) {
    return true;
  }
  boost::container::flat_set<int> const present(types.begin(), types.end());
  auto const pairs = std::views::cartesian_product(present, present);
  double const cutoff = std::ranges::max(
      pairs | std::views::transform([&](auto const &pair) {
        auto const [a, b] = pair;
        return tol.min_distance(a, b);
      }));

  MinimumImage const metric(lattice, periodicity);

  // Self-images: the nearest non-trivial image of any point is the same
  // lattice vector for every atom (infinite with no periodic axis).
  double const self_image =
      metric.distance(Vector3d::Zero(), Vector3d::Zero(), Images::nontrivial);
  if (std::ranges::any_of(present, [&](int type) {
        return self_image < tol.min_distance(type, type);
      })) {
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
      if (metric.distance(from, row(j), Images::all) <
          tol.min_distance(type_at(types, i), type_at(types, j))) {
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

bool distances_valid(Cell const &cell, DistanceTolerance const &tol) {
  return pairwise_distances_ok(cell.positions(), cell.types(),
                               cell.lattice().matrix(), cell.periodicity(),
                               tol);
}

} // namespace seitz::generate
