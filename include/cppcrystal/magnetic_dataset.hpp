#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <vector>

namespace cppcrystal::magnetic {

// Construction type of a magnetic space group (Barnighausen / BNS types I-IV).
enum class MagneticType { type_i = 1, type_ii = 2, type_iii = 3, type_iv = 4 };

} // namespace cppcrystal::magnetic

namespace cppcrystal {

// The result of a magnetic space-group determination (3D path). Carries the
// magnetic space-group identity, the magnetic operations of the input cell, the
// standardized cell (with site tensors), and the per-atom equivalence data.
struct MagneticDataset {
  UniNumber uni;
  magnetic::MagneticType msg_type = magnetic::MagneticType::type_i;
  // Hall number of the family (types I–III) or maximal (type IV) space group.
  HallNumber hall;
  SiteTensor tensor_rank = SiteTensor::collinear;

  // Magnetic operations of the input cell.
  MagneticOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;

  // Maps the input setting to the standardized one.
  Matrix3d transformation_matrix{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};

  // Standardized cell: idealized lattice, fractional positions, atom types, site
  // tensors (rotated into the standardized basis), and the rigid rotation that
  // orients the standardized lattice.
  Matrix3d std_lattice{Matrix3d::Identity()};
  Positions std_positions;
  std::vector<int> std_types;
  SiteTensors std_tensors;
  Matrix3d std_rotation_matrix{Matrix3d::Identity()};

  // Primitive lattice found during the magnetic symmetry search.
  Matrix3d primitive_lattice{Matrix3d::Identity()};
};

// Determine the magnetic space group of `cell` (positions + per-site magnetic
// tensors) and standardize it (3D path). The rank-1 tensors transform as axial
// or polar vectors according to `cell.kind()`; `tol.moment` (unset -> symprec)
// is the moment tolerance. Errors with e_magnetic_symmetry_search_failed /
// e_cell_standardization_failed on failure.
[[nodiscard]] Result<MagneticDataset>
get_magnetic_dataset(MagneticCell const &cell,
                     MagneticTolerance const &tol = {});

} // namespace cppcrystal
