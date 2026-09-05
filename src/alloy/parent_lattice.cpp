#include <seitz/alloy/parent_lattice.hpp>

#include "alloy/geometry.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>

#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <cstddef>
#include <format>
#include <functional>
#include <ranges>
#include <utility>
#include <vector>

namespace seitz::alloy {

using detail::integral;

Result<ParentLattice> ParentLattice::create(Cell cell,
                                            std::vector<Species> species,
                                            Tolerance tol) {
  if (species.empty()) {
    return leaf::new_error(
        e_message{"alloy: the parent lattice declares no sublattices"});
  }

  // "Ascending and without repeats" is one adjacent_find with a >= predicate:
  // it reports the first place the sequence fails to strictly increase, which
  // is exactly the disjunction of unsorted and duplicated.
  auto sublattices = std::views::enumerate(species);
  auto const malformed =
      std::ranges::find_if(sublattices, [](auto const &entry) {
        Species const &allowed = std::get<1>(entry);
        return allowed.empty() ||
               std::ranges::adjacent_find(allowed, std::ranges::greater_equal{}) !=
                   allowed.end();
      });
  if (malformed != sublattices.end()) {
    return leaf::new_error(e_message{std::format(
        "alloy: sublattice {} must admit at least one species, listed "
        "ascending and without repeats",
        std::get<0>(*malformed))});
  }

  auto atoms = std::views::enumerate(cell.types());
  auto const stray = std::ranges::find_if(atoms, [&](auto const &entry) {
    int const id = std::get<1>(entry);
    return id < 0 || std::cmp_greater_equal(id, species.size());
  });
  if (stray != atoms.end()) {
    return leaf::new_error(e_message{
        std::format("alloy: atom {} carries sublattice id {}, which has no "
                    "species entry",
                    std::get<0>(*stray), std::get<1>(*stray))});
  }

  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(cell, tol);
  BOOST_LEAF_AUTO(found, analyzer.cell_operations());
  // Copied out of the analyzer's memo: the analyzer is a local, and this is the
  // one thing the parent lattice keeps from it.
  Operations operations = found;
  return ParentLattice{std::move(cell), std::move(species),
                       std::move(operations)};
}

Result<ParentLattice> ParentLattice::from_sites(Lattice lattice,
                                                std::span<SiteSpec const> sites,
                                                Tolerance tol) {
  if (sites.empty()) {
    return leaf::new_error(e_message{"alloy: the parent lattice has no sites"});
  }

  // Sites admitting the same species set are one sublattice. Sorting each set
  // first makes "same chemistry" independent of the order the caller wrote it
  // in, which is what makes the derived ids -- and so the symmetry -- correct.
  boost::container::flat_map<Species, int> ids;
  std::vector<Species> species;
  std::vector<Vector3d> coordinates;
  Types types;
  coordinates.reserve(sites.size());
  types.reserve(sites.size());
  for (auto const &[position, allowed] : sites) {
    Species sorted = allowed;
    std::ranges::sort(sorted);
    auto const [entry, fresh] =
        ids.try_emplace(std::move(sorted), static_cast<int>(species.size()));
    if (fresh) {
      species.push_back(entry->first);
    }
    coordinates.push_back(position);
    types.push_back(entry->second);
  }
  return create(Cell{std::move(lattice), to_positions(coordinates),
                     std::move(types)},
                std::move(species), tol);
}

int ParentLattice::species_at(Index site) const noexcept {
  return species_of(cell_.type(site));
}

int ParentLattice::max_species() const noexcept {
  // fold_left from zero rather than ranges::max, which has an empty-range
  // precondition this would otherwise have to restate.
  return std::ranges::fold_left(
      species_ | std::views::transform([](Species const &allowed) {
        return static_cast<int>(allowed.size());
      }),
      0, [](int largest, int count) { return std::max(largest, count); });
}

std::vector<ClusterPoint> ParentLattice::anchors(Anchors which) const {
  // Cell::atoms() hands over each atom as {position, type}, so nothing here
  // indexes the cell and the species table in step.
  auto active = cell_.atoms() | std::views::filter([this](auto const &atom) {
                  return species_of(atom.second) > 1; // spectators carry none
                }) |
                std::views::transform([this](auto const &atom) {
                  return ClusterPoint{atom.first, species_of(atom.second), 0};
                });

  if (which == Anchors::all) {
    return std::vector<ClusterPoint>(active.begin(), active.end());
  }

  // The crystallographic orbits, computed from the operations already in hand
  // rather than from a full determination: a site is a duplicate when some
  // operation carries an already-kept site onto it modulo the cell.
  std::vector<ClusterPoint> representatives;
  for (ClusterPoint const &site : active) {
    bool const seen =
        std::ranges::any_of(representatives, [&](ClusterPoint const &kept) {
          return std::ranges::any_of(
              operations_, [&](SymmetryOperation const &op) {
                return integral(op.apply(kept.position) - site.position,
                                kClusterPrec);
              });
        });
    if (!seen) {
      representatives.push_back(site);
    }
  }
  return representatives;
}

Result<Cluster> ParentLattice::cluster_of(
    std::span<Vector3d const> points) const {
  Cluster::Points cluster;
  cluster.reserve(points.size());
  auto const atoms = cell_.atoms();
  for (auto const &[index, cartesian] : std::views::enumerate(points)) {
    Vector3d const fractional = lattice().to_fractional(cartesian);
    // Which site of the parent this point sits on, modulo the cell. The point
    // keeps its own coordinates: which periodic image it names is the cluster's
    // geometry, and only the species count comes from the site.
    auto const site = std::ranges::find_if(atoms, [&](auto const &atom) {
      return integral(fractional - atom.first, kClusterPrec);
    });
    if (site == std::ranges::end(atoms)) {
      return leaf::new_error(e_message{std::format(
          "alloy: maximal-cluster point {} lies on no site of the parent "
          "lattice",
          index)});
    }
    cluster.push_back(ClusterPoint{fractional, species_of((*site).second), 0});
  }
  return Cluster{std::move(cluster)};
}

} // namespace seitz::alloy
