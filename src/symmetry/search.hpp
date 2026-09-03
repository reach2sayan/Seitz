#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/tolerance.hpp>

#include <vector>

namespace cppcrystal::symmetry {

// The symmetry search over one cell at one tolerance: the lattice point group,
// the cell's own space-group operations, and the pure translations. The family
// is a compile-time parameter, so the layer path's 2D reduction and its halved
// operation cap are `if constexpr` branches rather than runtime tests.
//
// Non-owning: `cell` must outlive the search.
template <GroupFamily F> class SymmetrySearch {
public:
  SymmetrySearch(Cell const &cell, Tolerance const &tol) noexcept
      : cell_(cell), tol_(tol) {}

  // Point group of the lattice: the rotations (in the cell's own basis) that
  // map the Delaunay-reduced lattice metric onto itself. At most 48 (24 for a
  // layer). Errors with e_symmetry_operation_search_failed.
  [[nodiscard]] Result<PointSymmetry> lattice_symmetry() const;

  // All space-group operations of the cell exactly as given, including the
  // centering translations of a non-primitive cell. Errors with e_empty_cell
  // for a cell with no atoms, otherwise e_symmetry_operation_search_failed.
  [[nodiscard]] Result<Operations> operations() const;

  // Re-filter `operations` at a (usually tighter) tolerance: keep only those
  // whose rotation is still a lattice symmetry of the cell and whose
  // translation still maps it onto itself. Empty if the lattice symmetry
  // cannot be determined at that tolerance.
  [[nodiscard]] Operations reduce(Operations const &operations,
                                  Tolerance const &tol) const;

  // Every t for which (identity, t) maps the cell onto itself: the centering
  // translations, including the zero translation. These number
  // (cell size) / (primitive size).
  [[nodiscard]] std::vector<Vector3d> pure_translations() const;

private:
  Cell const &cell_;
  Tolerance tol_;
};

extern template class SymmetrySearch<GroupFamily::space>;
extern template class SymmetrySearch<GroupFamily::layer>;

} // namespace cppcrystal::symmetry
