#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/types.hpp>
#include <seitz/spacegroup_match.hpp>

#include <string_view>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::analysis {

// The per-atom result of a determination, one Site per input-cell atom.
struct Site {
  int wyckoff = 0;                // Wyckoff letter index, 0 = 'a'
  std::string_view site_symmetry; // tabulated site-symmetry symbol
  int equivalent_atom = 0;        // symmetrically-equivalent representative
  int orbit = 0;                  // crystallographic-orbit representative
  int primitive_atom = 0;         // this atom's primitive-cell atom
};

// The result of a space-group determination.
struct Dataset {
  HallNumber hall;
  Lattice bravais;
  Setting setting;

  // Space-group operations of the *input* cell (rotation in the cell basis,
  // fractional translation).
  Operations operations;
  std::vector<Site> sites;

  // Standardized conventional (bravais) cell -- idealized lattice, positions
  // and types with its own periodicity -- and each of its atoms' primitive-cell
  // atom.
  Cell standardized;
  std::vector<int> std_mapping_to_primitive;

  // Primitive cell's lattice, found during determination.
  Lattice primitive;
};

// The result of a magnetic space-group determination (3D path): the match
// (uni, type, hall, setting) plus the cell-level products.
struct MagneticDataset : MagneticMatch {
  MagneticOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;

  // Standardized magnetic cell, site tensors rotated into its basis.
  MagneticCell standardized;
  Lattice primitive;
};

} // namespace seitz::analysis

#pragma GCC visibility pop
