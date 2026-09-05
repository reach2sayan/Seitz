#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/types.hpp>

#include <string_view>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::magnetic {

// Construction type of a magnetic space group (Barnighausen / BNS types I-IV).
enum class MagneticType { type_i = 1, type_ii = 2, type_iii = 3, type_iv = 4 };

} // namespace seitz::magnetic

namespace seitz::analysis {

// input cell <--> standardized setting:
// (a) the change of basis + origin shift to align ops with the database,
// (b) rigid rotation to the idealized standardized lattice.
struct Setting {
  Matrix3d transformation{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};
  Matrix3d rigid_rotation{Matrix3d::Identity()};
};

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

  // Standardized conventional ("bravais") cell — idealized lattice, positions
  // and types, carrying its own periodicity — and the map from each of its
  // atoms to a primitive-cell atom.
  Cell standardized;
  std::vector<int> std_mapping_to_primitive;

  // Primitive cell's lattice, found during determination.
  Lattice primitive;
};

// The result of a magnetic space-group determination (3D path).
struct MagneticDataset {
  UniNumber uni;
  magnetic::MagneticType type = magnetic::MagneticType::type_i;
  // Family (types I-III) or maximal (type IV) space group.
  HallNumber hall;
  Setting setting;

  MagneticOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;

  // Standardized magnetic cell, site tensors rotated into its basis.
  MagneticCell standardized;
  Lattice primitive;
};

} // namespace seitz::analysis

#pragma GCC visibility pop
