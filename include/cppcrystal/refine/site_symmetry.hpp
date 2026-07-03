#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <string>
#include <vector>

namespace cppcrystal::refine {

// Per-atom exact-position / Wyckoff data for the conventional primitive cell.
struct ExactPositions {
  std::vector<Vector3d> positions;  // exact (symmetrized) fractional positions
  std::vector<int> wyckoffs;        // Wyckoff letter index (0 = a, 1 = b, ...)
  std::vector<int> equivalent_atoms;       // representative atom index per atom
  std::vector<std::string> site_symmetry_symbols;
};

// Exact positions + Wyckoff assignment for the atoms of `conv_prim` (positions
// expressed wrt the idealized conventional lattice). `conv_sym` are the
// conventional database operations of the Hall setting. Returns std::nullopt if
// the Wyckoff labels
// could not be resolved at any of the attempted tolerances.
[[nodiscard]] std::optional<ExactPositions>
exact_positions(Cell const &conv_prim, SymmetryOperations const &conv_sym,
                int num_pure_trans, int hall_number, double symprec);

} // namespace cppcrystal::refine
