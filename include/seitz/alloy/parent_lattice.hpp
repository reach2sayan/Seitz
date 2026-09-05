#pragma once

#include <seitz/alloy/cluster.hpp>
#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <boost/container/small_vector.hpp>

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::alloy {

// The atomic numbers one sublattice may hold, ascending. small_vector because
// every tractable expansion is binary or ternary; four inline covers
// quaternaries without touching the allocator in the loops that read it.
using Species = boost::container::small_vector<int, 4>;

// A site as a caller describes it: where it sits (fractional) and what may
// occupy it. The sublattice ids the Cell ends up carrying are derived from
// these, never supplied by hand.
using SiteSpec = std::pair<Vector3d, Species>;

// Which active sites seed cluster enumeration.
//
// `inequivalent` keeps one representative per crystallographic orbit: any
// cluster anchored on a site maps, under the group, onto one anchored on that
// site's representative, so the orbits generated are the same and the O(sites x
// candidates x orbits x |G|) dedup does strictly less work. `all` is the
// exhaustive alternative, kept as the escape hatch and as the reference the
// fast path is tested against.
enum class Anchors { inequivalent, all };

// The disordered parent of a cluster expansion, built on seitz::Cell.
//
// The Cell's types are SUBLATTICE IDS, not atomic numbers: two sites carry the
// same id exactly when they admit the same species set. That is precisely the
// equivalence the space-group search needs -- equal types means candidate
// images -- so the operations here are SymmetryAnalyzer's with nothing in
// between, and the per-site species-signature hash a cluster-expansion code
// normally needs disappears into Cell::types().
//
// cell_operations() rather than operations(): the enumeration wants every
// operation of the cell exactly as given, centering translations included, not
// the standardized dataset's. Using the latter would silently change every
// multiplicity this module reports.
//
// The operations are resolved once, in the factory, so the whole module has a
// single fallible entry point and every accessor below is infallible.
class ParentLattice {
public:
  // `cell.types()` are sublattice ids indexing `species`; `species[i]` lists
  // the atomic numbers allowed on sublattice i, ascending and without repeats.
  // A sublattice admitting one species is a spectator: it carries no cluster
  // point.
  [[nodiscard]] static Result<ParentLattice>
  create(Cell cell, std::vector<Species> species, Tolerance tol = {});

  // The constructor callers will actually use: one allowed-species list per
  // site. Sites with equal species sets are collapsed into one sublattice, ids
  // assigned in first-seen order, and create() does the rest.
  [[nodiscard]] static Result<ParentLattice>
  from_sites(Lattice lattice, std::span<SiteSpec const> sites,
             Tolerance tol = {});

  [[nodiscard]] Cell const &cell() const noexcept { return cell_; }
  [[nodiscard]] Lattice const &lattice() const noexcept {
    return cell_.lattice();
  }
  [[nodiscard]] Operations const &operations() const noexcept {
    return operations_;
  }
  // Species admitted by the sublattice of atom `site`; 1 means spectator.
  [[nodiscard]] int species_at(Index site) const noexcept;

  // Largest species count over all sublattices: the size of the site basis.
  [[nodiscard]] int max_species() const noexcept;

  // The sites that actually vary, as cluster points at lattice translation
  // zero with the constant point function. Every cluster enumeration starts
  // from these.
  [[nodiscard]] std::vector<ClusterPoint> anchors(Anchors which) const;

  // A maximal cluster given as CARTESIAN points, mapped into this lattice's
  // fractional frame with each point's species count filled in and every
  // function at zero -- which is exactly the shape Cvm::create expects. This is
  // the only place the module needs a site lookup, so it lives here rather than
  // leaving every caller to do it. Fails naming the point that lies on no site.
  [[nodiscard]] Result<Cluster> cluster_of(std::span<Vector3d const> points) const;

private:
  ParentLattice(Cell cell, std::vector<Species> species, Operations operations)
      : cell_{std::move(cell)}, species_{std::move(species)},
        operations_{std::move(operations)} {}

  // Species admitted by a sublattice id. The id-keyed form, so the enumeration
  // can walk Cell::atoms() -- which hands over the type directly -- instead of
  // indexing the cell and the table in step.
  [[nodiscard]] int species_of(int sublattice) const noexcept {
    return static_cast<int>(
        species_[static_cast<std::size_t>(sublattice)].size());
  }

  Cell cell_;
  // Keyed by the sublattice id in cell_.types(), which the factories make dense
  // from zero. A lookup table: nothing walks it in step with the sites.
  std::vector<Species> species_;
  Operations operations_;
};

} // namespace seitz::alloy

#pragma GCC visibility pop
