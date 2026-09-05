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

// The atomic numbers one sublattice may hold, ascending. Four inline covers
// every tractable (binary to quaternary) expansion allocation-free.
using Species = boost::container::small_vector<int, 4>;

// A site as a caller describes it: fractional position and admissible
// species. Sublattice ids are derived from these, never supplied.
using SiteSpec = std::pair<Vector3d, Species>;

// Which active sites seed cluster enumeration. `inequivalent` keeps one
// representative per crystallographic orbit -- every cluster maps onto one
// anchored there, so the orbit set is unchanged and the
// O(sites x candidates x orbits x |G|) dedup shrinks. `all` is the exhaustive
// reference the fast path is tested against.
enum class Anchors { inequivalent, all };

// The disordered parent of a cluster expansion, over seitz::Cell.
//
// Cell types are SUBLATTICE IDS, not atomic numbers: equal id iff equal species
// set, which is exactly the equivalence the space-group search wants, so the
// operations are SymmetryAnalyzer's unmediated.
//
// cell_operations(), not operations(): the enumeration needs every operation of
// the cell as given, centering translations included; the standardized
// dataset's would change every multiplicity reported here. Resolved once in the
// factory, so every accessor below is infallible.
class ParentLattice {
public:
  // `cell.types()` are sublattice ids indexing `species`; `species[i]` lists
  // sublattice i's atomic numbers, ascending and unique. A one-species
  // sublattice is a spectator and carries no cluster point.
  [[nodiscard]] static Result<ParentLattice>
  create(Cell cell, std::vector<Species> species, Tolerance tol = {});

  // One allowed-species list per site; equal species sets collapse into one
  // sublattice, ids in first-seen order.
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

  // The varying sites as cluster points at translation zero with the constant
  // point function; every enumeration starts from these.
  [[nodiscard]] std::vector<ClusterPoint> anchors(Anchors which) const;

  // A maximal cluster given as CARTESIAN points, mapped into the fractional
  // frame with species counts filled in and every function zero -- the shape
  // Cvm::create expects. Fails naming the point that lies on no site.
  [[nodiscard]] Result<Cluster> cluster_of(std::span<Vector3d const> points) const;

private:
  ParentLattice(Cell cell, std::vector<Species> species, Operations operations)
      : cell_{std::move(cell)}, species_{std::move(species)},
        operations_{std::move(operations)} {}

  // Species admitted by a sublattice id, so the enumeration can walk
  // Cell::atoms() directly instead of indexing cell and table in step.
  [[nodiscard]] int species_of(int sublattice) const noexcept {
    return static_cast<int>(
        species_[static_cast<std::size_t>(sublattice)].size());
  }

  Cell cell_;
  // Keyed by the sublattice id in cell_.types(), dense from zero.
  std::vector<Species> species_;
  Operations operations_;
};

} // namespace seitz::alloy

#pragma GCC visibility pop
