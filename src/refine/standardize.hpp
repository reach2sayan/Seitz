#pragma once

#include <seitz/analysis/dataset.hpp>
#include <seitz/core/cell.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#include <vector>

namespace seitz::refine {

// The standardized conventional ("bravais") cell together with the per-atom
// Wyckoff / equivalence data of the *input* cell.
struct Standardized {
  Cell bravais;                              // idealized conventional cell
  std::vector<analysis::Site> sites;                 // per input-cell atom
  std::vector<int> std_mapping_to_primitive; // per bravais atom
};

} // namespace seitz::refine
