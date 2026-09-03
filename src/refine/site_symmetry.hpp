#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <string>
#include <vector>

namespace cppcrystal::refine {

// Exact-position / Wyckoff data of one conventional-primitive atom.
struct ExactPosition {
  Vector3d position;        // exact (symmetrized) fractional position
  int wyckoff = 0;          // Wyckoff letter index (0 = a, 1 = b, ...)
  int equivalent_atom = 0;  // representative atom index
  std::string site_symmetry_symbol{};
};

// Per-atom exact-position / Wyckoff data for the conventional primitive cell.
using ExactPositions = std::vector<ExactPosition>;

// Exact positions + Wyckoff assignment for the atoms of `conv_prim` (positions
// expressed wrt the idealized conventional lattice). `conv_sym` are the
// conventional database operations of the Hall setting. Returns std::nullopt if
// the Wyckoff labels
// could not be resolved at any of the attempted tolerances.
[[nodiscard]] std::optional<ExactPositions>
exact_positions(Cell const &conv_prim, Operations const &conv_sym,
                int num_pure_trans, HallNumber hall, double symprec);

} // namespace cppcrystal::refine
