#include <seitz/alloy/clusters_pool.hpp>

#include "alloy/enumerate.hpp"
#include "alloy/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

namespace seitz::alloy {
namespace {

using detail::lattice_box;

// Every active point within Cartesian distance `radius` of `anchor`, the anchor
// excluded, over a box of integer translations sized from the metric.
//
// Brute force on purpose: EVERY neighbour inside the radius is wanted, several
// periodic images of one site included, so minimum_image_distance -- nearest
// image only -- cannot stand in. Sites outermost, which is the order the cutoff
// test sees.
[[nodiscard]] std::vector<ClusterPoint>
neighbours_within(Lattice const &lattice,
                  std::vector<ClusterPoint> const &active,
                  Vector3d const &anchor, double radius, double tol) {
  double const spread =
      radius * lattice.matrix().inverse().rowwise().norm().maxCoeff();
  int const reach = static_cast<int>(std::ceil(spread)) + 2;

  auto found =
      std::views::cartesian_product(active, lattice_box(reach)) |
      std::views::transform([](auto const &placement) {
        auto const &[site, translation] = placement;
        auto const [x, y, z] = translation;
        return ClusterPoint{site.position + Vector3d(x, y, z), site.species, 0};
      }) |
      std::views::filter([&](ClusterPoint const &candidate) {
        double const distance =
            lattice.to_cartesian(candidate.position - anchor).norm();
        return distance > tol && distance <= radius + tol;
      });
  return std::vector<ClusterPoint>(found.begin(), found.end());
}

// Record `candidate` as a new orbit unless an orbit equivalent to it is already
// present. The images are computed only for a genuinely new orbit, since
// generating them is the expensive half.
void add_orbit(Cluster candidate, std::vector<Orbit> &orbits,
               Operations const &ops, Lattice const &lattice, double tol) {
  if (find_equivalent(orbits, candidate, ops, tol, &Orbit::representative) !=
      orbits.end()) {
    return;
  }
  double const diameter = candidate.diameter(lattice);
  std::vector<Cluster> images = symmetry_images(candidate, ops, tol);
  orbits.push_back(Orbit{.representative = std::move(candidate),
                         .images = std::move(images),
                         .diameter = diameter});
}

// Every decoration of `geometry`, added to `orbits`. A site admitting k species
// carries k - 1 non-constant point functions, so that is the odometer's radix.
void add_decorations(Cluster const &geometry, std::vector<Orbit> &orbits,
                     Operations const &ops, Lattice const &lattice,
                     double tol) {
  std::vector<int> const radix = detail::radix_of(
      geometry, [](ClusterPoint const &p) { return p.species - 1; });
  for (std::span<int const> functions : detail::mixed_radix(radix)) {
    auto const decorated =
        std::views::zip(geometry, functions) |
        std::views::transform([](auto const &pair) {
          auto const &[point, function] = pair;
          return ClusterPoint{point.position, point.species, function};
        });
    add_orbit(Cluster{Cluster::Points(decorated.begin(), decorated.end())},
              orbits, ops, lattice, tol);
  }
}

} // namespace

Result<ClustersPool> ClustersPool::generate(ParentLattice const &parent,
                                            Options const &options) {
  Operations const &ops = parent.operations();
  Lattice const &lattice = parent.lattice();
  double const tol = options.tolerance;

  // Cluster enumeration is seeded from one site per crystallographic orbit, but
  // a cluster's OTHER points may sit on any active site, so the neighbour pool
  // is always the full set.
  std::vector<ClusterPoint> const seeds = parent.anchors(options.anchors);
  std::vector<ClusterPoint> const active = parent.anchors(Anchors::all);

  std::vector<Orbit> orbits;

  if (options.include_empty) {
    orbits.push_back(Orbit{.representative = {}, .images = {Cluster{}}});
  }

  if (options.include_points) {
    std::vector<Orbit> points;
    for (ClusterPoint const &seed : seeds) {
      add_decorations(Cluster{Cluster::Points{seed}}, points, ops, lattice,
                      tol);
    }
    std::ranges::move(points, std::back_inserter(orbits));
  }

  // Higher body orders, ascending (flat_map is sorted by key), each up to its
  // own diameter cutoff.
  for (auto const &[order, radius] : options.radii) {
    if (order < 2) {
      continue; // the empty and point clusters have their own switches
    }
    std::vector<Orbit> found;
    for (ClusterPoint const &seed : seeds) {
      std::vector<ClusterPoint> const nearby =
          neighbours_within(lattice, active, seed.position, radius, tol);
      for (std::span<std::size_t const> combination :
           detail::combinations(nearby.size(), order - 1)) {
        auto const chosen =
            combination | std::views::transform(
                              [&](std::size_t i) { return nearby[i]; });
        Cluster::Points geometry{seed};
        geometry.insert(geometry.end(), chosen.begin(), chosen.end());

        // The cutoff is on the cluster's own diameter, not on the distance from
        // the seed: a triangle of nearest neighbours can be wider than any of
        // its edges is long.
        Cluster const candidate{std::move(geometry)};
        if (candidate.diameter(lattice) > radius + tol) {
          continue;
        }
        add_decorations(candidate, found, ops, lattice, tol);
      }
    }
    // Stable, so orbits of equal diameter keep the order they were generated
    // in and the pool is reproducible run to run.
    std::ranges::stable_sort(found, std::ranges::less{}, &Orbit::diameter);
    std::ranges::move(found, std::back_inserter(orbits));
  }

  return ClustersPool{lattice, std::move(orbits)};
}

} // namespace seitz::alloy
