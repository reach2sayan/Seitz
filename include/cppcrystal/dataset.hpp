#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace cppcrystal {

// The result of a space-group determination (3D space-group path): the
// space-group identity, the operations of the input cell, the standardized
// conventional cell, and the per-atom Wyckoff / equivalence data.
struct Dataset {
  int spacegroup_number = 0; // international / layer-group number
  HallNumber hall;           // the matched Hall setting

  std::string_view international_symbol; // Hermann-Mauguin (short)
  std::string_view hall_symbol;
  std::string_view choice;

  int pointgroup_number = 0; // 1..32
  std::string_view pointgroup_symbol;

  // Maps the input cell onto the conventional/bravais setting:
  //   bravais_lattice = cell.lattice . transformation_matrix^-1
  // and the origin shift aligning the operations with the Hall database.
  Matrix3d bravais_lattice{Matrix3d::Identity()};
  Matrix3d transformation_matrix{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};

  // Space-group operations of the *input* cell (rotation in the cell basis,
  // fractional translation).
  Operations operations;

  // Per input-cell atom: Wyckoff letter index (0 = a), site-symmetry symbol,
  // symmetrically-equivalent representative atom, and crystallographic orbit
  // representative (these differ only for symmetry-broken supercells).
  std::vector<int> wyckoffs;
  std::vector<std::string> site_symmetry_symbols;
  std::vector<int> equivalent_atoms;
  std::vector<int> crystallographic_orbits;

  // Standardized conventional ("bravais") cell: idealized lattice, fractional
  // positions and atom types, the rigid rotation that orients the standardized
  // basis, and the map from each standardized atom to its primitive atom.
  Matrix3d std_lattice{Matrix3d::Identity()};
  Positions std_positions;
  std::vector<int> std_types;
  Matrix3d std_rotation_matrix{Matrix3d::Identity()};
  std::vector<int> std_mapping_to_primitive;

  // Primitive cell found during determination and the map from each input-cell
  // atom to its primitive atom.
  Matrix3d primitive_lattice{Matrix3d::Identity()};
  std::vector<int> mapping_to_primitive;

  // The aperiodic axis (0/1/2) for a layer-group result, or std::nullopt for a
  // 3D space-group result. In the standardized cell this is always axis 2 (c).
  std::optional<int> aperiodic_axis;
};

// Determine the space group of `cell` and standardize it (3D space-group path):
// find_primitive -> search_spacegroup -> refinement (standardized cell + Wyckoff
// data). An unset `setting` searches every Hall setting of the cell's family.
//
// Errors with e_spacegroup_search_failed / e_cell_standardization_failed when
// determination fails at every attempted tolerance.
//
// For the object-oriented entry point that owns the cell and memoizes this and
// every derived query, prefer cppcrystal::analysis::SymmetryAnalyzer::dataset().
[[nodiscard]] Result<Dataset>
get_dataset(Cell const &cell, Tolerance const &tol = {},
            std::optional<HallNumber> setting = std::nullopt);

} // namespace cppcrystal
