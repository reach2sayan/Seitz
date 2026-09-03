#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>

// For the object-oriented entry points over a single cell, prefer
// cppcrystal::analysis::SymmetryAnalyzer::lattice_symmetry() and ::cell_operations(),
// which own the cell and memoize these results.
namespace cppcrystal::symmetry {

// Point group of the lattice: the rotations (in the cell's own basis) that map
// the Delaunay-reduced lattice metric onto itself. At most 48. Errors with
// e_symmetry_operation_search_failed.
[[nodiscard]] Result<PointSymmetry> lattice_symmetry(Cell const &cell,
                                                     Tolerance const &tol);

// All space-group operations of the cell exactly as given (including the
// centering/pure translations of a non-primitive cell). Errors with
// e_empty_cell for a cell with no atoms, otherwise
// e_symmetry_operation_search_failed.
[[nodiscard]] Result<SymmetryOperations> find_symmetry(Cell const &cell,
                                                       Tolerance const &tol);

// Re-filter `operations` at a (usually tighter) tolerance: keep only those
// whose rotation is still a lattice symmetry of `cell` at `symprec` and whose
// translation still maps the cell onto itself within `symprec`. Returns an empty
// set if the lattice symmetry cannot be determined at this tolerance.
[[nodiscard]] SymmetryOperations
reduce_symmetry(Cell const &cell, SymmetryOperations const &operations,
                Tolerance const &tol);

// Pure translations: every t for which (identity, t) maps the cell onto itself
// (the centering translations, including the zero translation). These number
// (cell size) / (primitive size).
[[nodiscard]] std::vector<Vector3d> pure_translations(Cell const &cell,
                                                      double symprec);

} // namespace cppcrystal::symmetry
