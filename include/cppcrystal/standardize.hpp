#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/tolerance.hpp>

namespace cppcrystal {

// Which setting the standardized cell is expressed in.
enum class CellSetting { conventional, primitive };

// Whether the standardized cell carries the metric-idealized (symmetrized)
// lattice, or keeps the input's real — possibly distorted — geometry and is
// only transformed into the standardized basis/centering.
enum class Idealize { yes, no };

// Standardize `cell` (3D space-group path). The four combinations:
//   {conventional, yes} the dataset's std_* cell
//   {primitive,    yes} the idealized primitive cell
//   {conventional, no } conventional basis, input geometry preserved
//   {primitive,    no } primitive basis, input geometry preserved
// Errors with e_spacegroup_search_failed / e_cell_standardization_failed when
// determination fails (same conditions as get_dataset).
//
// For the object-oriented entry point, prefer
// cppcrystal::analysis::SymmetryAnalyzer::standardized_cell(...).
[[nodiscard]] Result<Cell>
standardize_cell(Cell const &cell,
                 CellSetting setting = CellSetting::conventional,
                 Idealize idealize = Idealize::yes, Tolerance const &tol = {});

} // namespace cppcrystal
