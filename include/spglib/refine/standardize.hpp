#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>
#include <spglib/spacegroup/spacegroup.hpp>

#include <optional>
#include <string>
#include <vector>

namespace spglib::refine {

// The standardized conventional ("bravais") cell together with the per-atom
// Wyckoff / equivalence data of the *input* cell. Port of the get_Wyckoff_
// positions assembly in refinement.c (3D path).
struct Standardized {
  Cell bravais;                              // idealized conventional cell
  std::vector<int> wyckoffs;                 // per input-cell atom (0 = a, ...)
  std::vector<std::string> site_symmetry_symbols; // per input-cell atom
  std::vector<int> equivalent_atoms;         // per input-cell atom
  std::vector<int> crystallographic_orbits;  // per input-cell atom
  std::vector<int> std_mapping_to_primitive; // per bravais atom
};

// Build the standardized cell and per-atom Wyckoff data. `sg` must already be
// adjusted by find_similar_bravais_lattice. `cell_operations` are the symmetry
// operations of the input cell (used only for the supercell broken-symmetry
// check / equivalent-atom search). `mapping_table[i]` maps input-cell atom i to
// its primitive-cell atom. Returns std::nullopt if the Wyckoff resolution
// failed. Port of refinement.c get_Wyckoff_positions (+ helpers).
[[nodiscard]] std::optional<Standardized>
get_wyckoff_positions(spacegroup::Spacegroup const &sg, Cell const &primitive,
                      Cell const &cell,
                      SymmetryOperations const &cell_operations,
                      std::vector<int> const &mapping_table, double symprec);

} // namespace spglib::refine
