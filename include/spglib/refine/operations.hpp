#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/spacegroup/spacegroup.hpp>

#include <optional>

namespace spglib::refine {

// Exact, database-derived space-group operations of the input cell: the
// conventional database operations of the Hall setting, shifted by the origin
// shift, transformed to the primitive setting, then recovered in the input
// cell (replicated across its pure translations). Port of refinement.c
// get_refined_symmetry_operations (+ recover_symmetry_in_original_cell).
//
// `sg` must already be adjusted by find_similar_bravais_lattice. Returns
// std::nullopt when the recovered pure-translation count is inconsistent with
// the input cell's multiplicity (the caller retries at a tighter tolerance).
// 3D space-group path only.
[[nodiscard]] std::optional<SymmetryOperations>
refined_operations(spacegroup::Spacegroup const &sg, Cell const &primitive,
                   Cell const &cell, double symprec);

} // namespace spglib::refine
