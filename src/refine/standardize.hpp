#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <string_view>
#include <vector>

namespace cppcrystal::refine {

// The standardized conventional ("bravais") cell together with the per-atom
// Wyckoff / equivalence data of the *input* cell.
struct Standardized {
  Cell bravais;                              // idealized conventional cell
  std::vector<int> wyckoffs;                 // per input-cell atom (0 = a, ...)
  // Views into the static site-symmetry table, one per input-cell atom.
  std::vector<std::string_view> site_symmetry_symbols;
  std::vector<int> equivalent_atoms;         // per input-cell atom
  std::vector<int> crystallographic_orbits;  // per input-cell atom
  std::vector<int> std_mapping_to_primitive; // per bravais atom
};

} // namespace cppcrystal::refine
