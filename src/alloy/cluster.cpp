#include <cppcrystal/alloy/cluster.hpp>

#include "alloy/geometry.hpp"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <ranges>
#include <utility>
#include <vector>

namespace cppcrystal::alloy {

using detail::coincident;
using detail::integral;

Cluster Cluster::transformed(SymmetryOperation const &op) const {
  auto const moved =
      points_ | std::views::transform([&op](ClusterPoint point) {
        point.position = op.apply(point.position);
        return point;
      });
  return Cluster{Points(moved.begin(), moved.end())};
}

double Cluster::diameter(Lattice const &lattice) const {
  // Max over unordered pairs. Max is exact, so unlike a sum the fold order
  // does not affect the result bit for bit.
  auto const indices = std::views::iota(0, size());
  // Not const: filter_view caches its first element, so its begin() is
  // non-const and a const filter_view is not a range at all.
  auto pairs = std::views::cartesian_product(indices, indices) |
               std::views::filter([](auto const &pair) {
                 return std::get<0>(pair) < std::get<1>(pair);
               });
  return std::ranges::fold_left(
      pairs, 0.0, [&](double longest, auto const &pair) {
        auto const [first, second] = pair;
        Vector3d const delta =
            (*this)[first].position - (*this)[second].position;
        return std::max(longest, lattice.to_cartesian(delta).norm());
      });
}

bool Cluster::congruent(Cluster const &other, double tol) const {
  if (size() != other.size()) {
    return false;
  }
  if (empty()) {
    return true;
  }
  // Any point of `other` congruent to our first point modulo the cell fixes a
  // candidate whole-cluster shift; the cluster matches if that shift carries
  // every one of our points onto a point of `other` with the same function.
  Vector3d const first = points_.front().position;
  return std::ranges::any_of(other, [&](ClusterPoint const &anchor) {
    if (!integral(first - anchor.position, tol)) {
      return false;
    }
    Vector3d const shift = anchor.position - first;
    return std::ranges::all_of(points_, [&](ClusterPoint const &point) {
      Vector3d const target = point.position + shift;
      return std::ranges::any_of(other, [&](ClusterPoint const &candidate) {
        return candidate.function == point.function &&
               (candidate.position - target).norm() < tol;
      });
    });
  });
}

Cell Cluster::as_cell(Lattice const &lattice) const {
  auto const coordinates =
      points_ | std::views::transform(&ClusterPoint::position);
  std::vector<Vector3d> const rows(coordinates.begin(), coordinates.end());
  return Cell{lattice, to_positions(rows),
              Types{std::from_range,
                    points_ | std::views::transform(&ClusterPoint::function)}};
}

std::vector<Cluster> symmetry_images(Cluster const &c, Operations const &ops,
                                     double tol) {
  // Sequential by nature: each image is kept only if it is distinct from the
  // ones already kept, so the result is an accumulation and not a projection.
  std::vector<Cluster> images;
  for (SymmetryOperation const &op : ops) {
    Cluster image = c.transformed(op);
    if (std::ranges::none_of(images, [&](Cluster const &seen) {
          return seen.congruent(image, tol);
        })) {
      images.push_back(std::move(image));
    }
  }
  return images;
}

Operations site_symmetry(Cluster const &c, Operations const &ops, double tol) {
  std::vector<SymmetryOperation> stabilizer;
  for (SymmetryOperation const &op : ops) {
    Cluster const image = c.transformed(op);
    if (!image.congruent(c, tol)) {
      continue;
    }
    // All points share one lattice shift, so their mean residual is that shift;
    // subtracting it leaves an operation that permutes the cluster's own points
    // with no net translation, which is what point_permutation needs.
    Vector3d const residual =
        std::ranges::fold_left(
            std::views::zip(c, image), Vector3d{Vector3d::Zero()},
            [](Vector3d sum, auto const &pair) {
              auto const &[original, moved] = pair;
              return Vector3d(sum + original.position - moved.position);
            }) /
        static_cast<double>(c.size());
    stabilizer.push_back(
        SymmetryOperation{op.rotation, op.translation + residual});
  }
  return Operations{std::move(stabilizer)};
}

std::vector<int> point_permutation(Cluster const &c,
                                   SymmetryOperation const &op, double tol) {
  Cluster const image = c.transformed(op);
  auto const targets =
      std::views::enumerate(image) |
      std::views::transform([&](auto const &pair) {
        auto const &[index, moved] = pair;
        auto const at = coincident(c, moved.position, tol);
        // site_symmetry guarantees a hit; a miss would mean the op does not
        // stabilize the cluster, so map the point to itself rather than
        // invent a target.
        return at == c.end()
                   ? static_cast<int>(index)
                   : static_cast<int>(std::ranges::distance(c.begin(), at));
      });
  return std::vector<int>{std::from_range, targets};
}

} // namespace cppcrystal::alloy
