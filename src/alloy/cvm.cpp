#include <seitz/alloy/cvm.hpp>

#include "alloy/enumerate.hpp"
#include "alloy/geometry.hpp"

#include <algorithm>
#include <cstddef>
#include <format>
#include <functional>
#include <limits>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

namespace seitz::alloy {
namespace {

using detail::coincident;
using detail::integral;
using detail::lattice_box;
using detail::mixed_radix;
using detail::radix_of;

// Every subset of the maximal cluster's ACTIVE points, deduplicated by the full
// space group and inserted so the list stays ordered by (point count,
// diameter). Inserting rather than collecting-and-sorting is what makes a tie
// keep the FIRST-generated representative, which the ordering of everything
// downstream -- and the Kikuchi-Barker recursion in particular -- depends on.
//
// The list accumulates across all the maximal clusters, so a subcluster two of
// them share is recorded once.
void collect_subclusters(std::vector<Cluster> &subclusters,
                         Cluster const &maximal, Operations const &ops,
                         Lattice const &lattice, double tol) {
  // A spectator point is never in a subcluster, so its digit has one value.
  std::vector<int> const radix = radix_of(
      maximal, [](ClusterPoint const &p) { return p.species > 1 ? 2 : 1; });

  for (std::span<int const> selection : mixed_radix(radix)) {
    auto taken = std::views::zip(selection, maximal) |
                 std::views::filter([](auto const &pair) {
                   return std::get<0>(pair) == 1;
                 }) |
                 std::views::transform([](auto const &pair) {
                   ClusterPoint const &point = std::get<1>(pair);
                   return ClusterPoint{point.position, point.species, 0};
                 });
    Cluster candidate{Cluster::Points(taken.begin(), taken.end())};
    if (find_equivalent(subclusters, candidate, ops, tol) !=
        subclusters.end()) {
      continue;
    }
    double const length = candidate.diameter(lattice);
    auto const slot =
        std::ranges::find_if(subclusters, [&](Cluster const &present) {
          if (present.size() != candidate.size()) {
            return present.size() > candidate.size();
          }
          return present.diameter(lattice) > length;
        });
    subclusters.insert(slot, std::move(candidate));
  }
}

// Every point-function decoration of every subcluster, deduplicated within its
// own subcluster by the full space group. A binary alloy yields one decoration
// per subcluster; a ternary point yields two.
[[nodiscard]] std::vector<ClusterFunction>
build_functions(std::span<Cluster const> subclusters, Operations const &ops,
                double tol) {
  std::vector<ClusterFunction> functions;
  for (Cluster const &sub : subclusters) {
    std::vector<int> const radix =
        radix_of(sub, [](ClusterPoint const &p) { return p.species - 1; });

    std::vector<Cluster> distinct;
    for (std::span<int const> assignment : mixed_radix(radix)) {
      auto const decorated =
          std::views::zip(sub, assignment) |
          std::views::transform([](auto const &pair) {
            auto const &[point, function] = pair;
            return ClusterPoint{point.position, point.species, function};
          });
      Cluster candidate{Cluster::Points(decorated.begin(), decorated.end())};
      if (find_equivalent(distinct, candidate, ops, tol) == distinct.end()) {
        distinct.push_back(std::move(candidate));
      }
    }
    auto decorated =
        distinct | std::views::transform([&](Cluster const &candidate) {
          return ClusterFunction{.cluster = candidate,
                                 .images =
                                     symmetry_images(candidate, ops, tol)};
        });
    functions.insert(functions.end(), decorated.begin(), decorated.end());
  }
  return functions;
}

// Every species assignment on `sub`, grouped into orbits of the subcluster's
// own site symmetry.
//
// The grouping compares occupations through PERMUTATIONS of the point indices
// rather than by rebuilding and geometrically matching a decorated cluster for
// every candidate. Each stabilizer element permutes the subcluster's own points
// -- that is what the translation correction in site_symmetry() buys -- so the
// permutations are computed once and equivalence becomes ranges::equal against
// the permuted occupation.
[[nodiscard]] std::vector<Configuration>
build_configurations(Cluster const &sub, Operations const &ops, double tol) {
  Operations const stabilizer = site_symmetry(sub, ops, tol);
  auto const permuted_by =
      stabilizer | std::views::transform([&](SymmetryOperation const &op) {
        return point_permutation(sub, op, tol);
      });
  std::vector<std::vector<int>> const permutations{std::from_range,
                                                   permuted_by};

  std::vector<int> const radix =
      radix_of(sub, [](ClusterPoint const &p) { return p.species; });

  std::vector<Configuration> configurations;
  for (std::span<int const> occupation : mixed_radix(radix)) {
    auto const match = std::ranges::find_if(
        configurations, [&](Configuration const &present) {
          return std::ranges::any_of(
              permutations, [&](std::vector<int> const &permutation) {
                auto const relabelled =
                    permutation | std::views::transform([&](int target) {
                      return present
                          .occupation[static_cast<std::size_t>(target)];
                    });
                return std::ranges::equal(occupation, relabelled);
              });
        });
    if (match == configurations.end()) {
      configurations.push_back(
          Configuration{.occupation = {occupation.begin(), occupation.end()},
                        .multiplicity = 1});
    } else {
      ++match->multiplicity;
    }
  }
  return configurations;
}

// The componentwise range of a set of candidate lattice translations.
struct TranslationBox {
  Vector3d low{Vector3d::Constant(std::numeric_limits<double>::max())};
  Vector3d high{Vector3d::Constant(std::numeric_limits<double>::lowest())};
};

// The correlation of one species assignment on `sites` under a set of cluster
// images: the sum, over every lattice placement of every image that lands
// entirely on the assignment's points, of the product of the point functions.
//
// A raw sum, not an average -- the normalization is the v-matrix's prefactor.
[[nodiscard]] double
configuration_correlation(Cluster const &sites,
                          std::span<int const> occupation,
                          std::span<Cluster const> images,
                          SiteBasis const &basis, double tol) {
  if (images.empty() || images.front().empty()) {
    return 1.0; // the empty cluster function correlates to one by definition
  }
  if (sites.empty()) {
    // The empty subcluster: no point for an image point to land on, so no
    // placement contributes. Its one configuration then has probability
    // V . xi = 1 * xi_empty = 1, which is the identity the whole expansion is
    // normalized against. (The reference implementation misses this case and
    // sizes its translation box off an empty range.)
    return 0.0;
  }

  // The translations worth trying are those aligning some image point onto some
  // occupied point. In the parent's fractional frame a lattice translation is
  // an integer triple, so the box is the componentwise range of the rounded
  // differences over every (image, occupied point, image point) triple.
  // Not const: the joined inner views are prvalues, so join_view keeps the
  // current one in a member and cannot iterate through a const reference.
  auto offsets =
      images | std::views::transform([&sites](Cluster const &image) {
        return std::views::cartesian_product(sites, image) |
               std::views::transform([](auto const &pair) {
                 auto const &[occupied, point] = pair;
                 return Vector3d((occupied.position - point.position)
                                     .array()
                                     .round()
                                     .matrix());
               });
      }) |
      std::views::join;
  auto const box = std::ranges::fold_left(
      offsets, TranslationBox{}, [](TranslationBox span, Vector3d const &offset) {
        return TranslationBox{Vector3d(span.low.cwiseMin(offset)),
                              Vector3d(span.high.cwiseMax(offset))};
      });

  // Occupations are looked up through the zip, so the assignment and the
  // geometry it decorates are never indexed in step.
  auto occupied = std::views::zip(sites, occupation);

  // Every image against every candidate translation, images outermost -- the
  // order the sum accumulates in, and so part of the result.
  Vector3i const first = box.low.cast<int>();
  Vector3i const last = box.high.cast<int>();
  auto placements =
      std::views::cartesian_product(images, lattice_box(first, last));

  return std::ranges::fold_left(
      placements, 0.0, [&](double total, auto const &placement) {
        auto const &[image, translation] = placement;
        auto const [x, y, z] = translation;
        Vector3d const shift(x, y, z);

        double product = 1.0;
        // The product short-circuits: an image point landing on no occupied
        // point drops the whole placement, so `product` is only ever
        // accumulated for placements that fully land.
        bool const landed =
            std::ranges::all_of(image, [&](ClusterPoint const &point) {
              Vector3d const target = point.position + shift;
              auto const at =
                  std::ranges::find_if(occupied, [&](auto const &entry) {
                    return (std::get<0>(entry).position - target).norm() < tol;
                  });
              if (at == occupied.end()) {
                return false;
              }
              product *= basis[point.species, point.function, std::get<1>(*at)];
              return true;
            });
        return landed ? total + product : total;
      });
}

// The v-matrix of one subcluster: rows are its configurations, columns are ALL
// cluster functions (columns for functions that cannot overlap this subcluster
// come out zero, which is what makes rho = V . xi a plain matrix product over
// the whole correlation vector).
[[nodiscard]] MatrixXd
build_vmatrix(Cluster const &sub,
              std::span<Configuration const> configurations,
              std::span<ClusterFunction const> functions,
              SiteBasis const &basis, double tol) {
  MatrixXd vmatrix(static_cast<Index>(configurations.size()),
                   static_cast<Index>(functions.size()));
  for (auto const &[column, function] : std::views::enumerate(functions)) {
    // The normalization that turns correlations into probabilities: a uniform
    // prior (the species count) on each site the function does not touch, and
    // the squared norm of the point-function row on each site it does.
    double const prefactor = std::ranges::fold_left(
        sub, 1.0, [&](double running, ClusterPoint const &point) {
          auto const at = coincident(function.cluster, point.position, tol);
          return running * (at == function.cluster.end()
                                ? static_cast<double>(point.species)
                                : basis.sum_of_squares(at->species,
                                                       at->function));
        });
    for (auto const &[row, configuration] :
         std::views::enumerate(configurations)) {
      vmatrix(row, column) =
          configuration_correlation(sub, configuration.occupation,
                                    function.images, basis, tol) /
          prefactor;
    }
  }
  return vmatrix;
}

// One overlap of a subcluster with a site: the re-anchored cluster and the
// index of the subcluster it came from.
struct Overlap {
  Cluster cluster;
  std::size_t subcluster = 0;
};

// Every symmetry image of every non-empty subcluster that, re-anchored on
// `site`, passes through it. The overlapping point is moved to index 0 and the
// rest follow, carried by the same lattice shift.
[[nodiscard]] std::vector<Overlap>
overlaps_at_site(Vector3d const &site, std::span<Cluster const> subclusters,
                 Operations const &ops, double tol) {
  std::vector<Overlap> overlaps;
  for (auto const &[index, sub] : std::views::enumerate(subclusters)) {
    if (sub.empty()) {
      continue;
    }
    for (Cluster const &image : symmetry_images(sub, ops, tol)) {
      for (auto const &[centre, anchor] : std::views::enumerate(image)) {
        Vector3d const shift = site - anchor.position;
        // A plain integrality test on the fractional shift. The reference
        // implementation multiplies by the inverse cell first, which is a
        // second change of basis on coordinates that are already fractional; it
        // agrees with this only because the fcc case it is tuned on happens to
        // have an integer inverse.
        if (!integral(shift, tol)) {
          continue;
        }
        auto carried = std::views::enumerate(image) |
                       std::views::filter([centre](auto const &pair) {
                         return std::get<0>(pair) != centre;
                       }) |
                       std::views::transform([&shift](auto const &pair) {
                         ClusterPoint point = std::get<1>(pair);
                         point.position += shift;
                         return point;
                       });
        Cluster::Points moved{
            ClusterPoint{site, anchor.species, anchor.function}};
        moved.insert(moved.end(), carried.begin(), carried.end());
        overlaps.push_back(
            Overlap{.cluster = Cluster{std::move(moved)},
                    .subcluster = static_cast<std::size_t>(index)});
      }
    }
  }
  return overlaps;
}

// The Kikuchi-Barker coefficients: the Moebius inversion over the subcluster
// inclusion lattice,
//
//   k_c = 1 - sum over strictly larger clusters c' containing c of k_c',
//
// with the empty cluster at zero. Walked largest index to smallest, which is
// valid only because the subcluster list is ordered by point count -- every
// term on the right is already final by the time it is read.
[[nodiscard]] std::vector<double>
build_kikuchi_barker(std::span<Cluster const> subclusters,
                     Operations const &ops, double tol) {
  std::vector<double> coefficients(subclusters.size(), 0.0);
  for (std::size_t const current :
       std::views::iota(std::size_t{1}, subclusters.size()) |
           std::views::reverse) {
    Cluster const &sub = subclusters[current];

    auto containing =
        overlaps_at_site(sub[0].position, subclusters, ops, tol) |
        std::views::filter([&](Overlap const &overlap) {
          return overlap.cluster.size() > sub.size() && // strictly larger
                 std::ranges::all_of(sub, [&](ClusterPoint const &point) {
                   return coincident(overlap.cluster, point.position, tol) !=
                          overlap.cluster.end();
                 });
        });
    coefficients[current] = std::ranges::fold_left(
        containing, 1.0, [&](double running, Overlap const &overlap) {
          return running - coefficients[overlap.subcluster];
        });
  }
  return coefficients;
}

} // namespace

Result<Cvm> Cvm::create(ParentLattice const &parent,
                        std::span<Cluster const> maximal,
                        SiteBasis const &basis, double tolerance) {
  if (maximal.empty()) {
    return leaf::new_error(
        e_message{"alloy: the CVM needs at least one maximal cluster"});
  }
  if (basis.max_species() < parent.max_species()) {
    return leaf::new_error(e_message{std::format(
        "alloy: the site basis covers up to {} species but the parent lattice "
        "has a sublattice admitting {}",
        basis.max_species(), parent.max_species())});
  }

  Operations const &ops = parent.operations();
  Lattice const &lattice = parent.lattice();

  std::vector<Cluster> subclusters;
  for (Cluster const &cluster : maximal) {
    collect_subclusters(subclusters, cluster, ops, lattice, tolerance);
  }

  std::vector<ClusterFunction> functions =
      build_functions(subclusters, ops, tolerance);
  std::vector<double> const coefficients =
      build_kikuchi_barker(subclusters, ops, tolerance);

  std::vector<CvmCluster> clusters;
  clusters.reserve(subclusters.size());
  for (auto &&[sub, coefficient] : std::views::zip(subclusters, coefficients)) {
    std::vector<Configuration> configurations =
        build_configurations(sub, ops, tolerance);
    MatrixXd vmatrix =
        build_vmatrix(sub, configurations, functions, basis, tolerance);
    clusters.push_back(CvmCluster{.sites = sub,
                                  .configurations = std::move(configurations),
                                  .vmatrix = std::move(vmatrix),
                                  .kikuchi_barker = coefficient,
                                  .diameter = sub.diameter(lattice)});
  }

  return Cvm{lattice, std::move(clusters), std::move(functions)};
}

} // namespace seitz::alloy
