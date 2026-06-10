#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/point_group.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>

namespace spglib::symmetry {

// Point group of the lattice: the rotations (in the cell's own basis) that map
// the Delaunay-reduced lattice metric onto itself. At most 48. Port of
// symmetry.c get_lattice_symmetry. Errors with
// e_symmetry_operation_search_failed.
[[nodiscard]] Result<PointSymmetry>
lattice_symmetry(Cell const &cell, double symprec,
                 AngleTolerance angle_tolerance = std::nullopt);

// All space-group operations of the cell exactly as given (including the
// centering/pure translations of a non-primitive cell). Port of
// symmetry.c sym_get_operation / get_space_group_operations. Errors with
// e_empty_cell for a cell with no atoms, otherwise
// e_symmetry_operation_search_failed.
[[nodiscard]] Result<SymmetryOperations>
find_symmetry(Cell const &cell, double symprec,
              AngleTolerance angle_tolerance = std::nullopt);

// Re-filter `operations` at a (usually tighter) tolerance: keep only those
// whose rotation is still a lattice symmetry of `cell` at `symprec` and whose
// translation still maps the cell onto itself within `symprec`. Port of
// symmetry.c sym_reduce_operation. Returns an empty set if the lattice symmetry
// cannot be determined at this tolerance.
[[nodiscard]] SymmetryOperations
reduce_symmetry(Cell const &cell, SymmetryOperations const &operations,
                double symprec, AngleTolerance angle_tolerance = std::nullopt);

// Pure translations: every t for which (identity, t) maps the cell onto itself
// (the centering translations, including the zero translation). These number
// (cell size) / (primitive size). Port of sym_get_pure_translation.
[[nodiscard]] std::vector<Vector3d> pure_translations(Cell const &cell,
                                                      double symprec);

} // namespace spglib::symmetry
