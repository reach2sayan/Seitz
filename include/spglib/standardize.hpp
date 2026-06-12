#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>

namespace spglib {

// Flags controlling cell standardization (spglib spg_standardize_cell).
struct StandardizeOptions {
  // Return the standardized primitive cell instead of the conventional one.
  bool to_primitive = false;
  // Keep the input cell's actual (possibly distorted) lattice rather than the
  // metric-idealized bravais lattice: the input is only transformed into the
  // standardized basis/centering, not symmetrized/idealized.
  bool no_idealize = false;
};

// Standardize `cell` to its conventional or primitive setting. Port of spglib's
// spg_standardize_cell (3D space-group path). The four flag combinations:
//   {false, false} conventional, idealized   (== the dataset std_* cell)
//   {true,  false} primitive,    idealized
//   {false, true } conventional, input geometry preserved (no idealization)
//   {true,  true } primitive,    input geometry preserved
// Errors with e_spacegroup_search_failed / e_cell_standardization_failed when
// determination fails (same conditions as get_dataset).
//
// For the object-oriented entry point, prefer
// spglib::analysis::SymmetryAnalyzer::standardized_cell(options).
[[nodiscard]] Result<Cell>
standardize_cell(Cell const &cell, StandardizeOptions options = {},
                 double symprec = kDefaultSymprec,
                 AngleTolerance angle_tolerance = std::nullopt);

// Symmetrize `cell` into its idealized conventional setting — spglib's
// spg_refine_cell, exactly standardize_cell with the default options.
[[nodiscard]] inline Result<Cell>
refine_cell(Cell const &cell, double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt) {
  return standardize_cell(cell, {}, symprec, angle_tolerance);
}

// Reduce `cell` to its standardized (idealized) primitive setting — spglib's
// spg_find_primitive, exactly standardize_cell with to_primitive set. Note this
// is distinct from symmetry::find_primitive (= prm_get_primitive), which returns
// the Primitive struct rather than the standardized cell.
[[nodiscard]] inline Result<Cell>
find_primitive(Cell const &cell, double symprec = kDefaultSymprec,
               AngleTolerance angle_tolerance = std::nullopt) {
  return standardize_cell(cell, {.to_primitive = true}, symprec,
                          angle_tolerance);
}

} // namespace spglib
