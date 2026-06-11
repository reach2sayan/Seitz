#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/magnetic/magnetic_spacegroup.hpp> // magnetic::MagneticType

#include <optional>
#include <vector>

namespace spglib {

// The result of a magnetic space-group determination — the modern-C++ value
// analogue of spglib's SpglibMagneticDataset (3D path). Carries the magnetic
// space-group identity, the magnetic operations of the input cell, the
// standardized cell (with site tensors), and the per-atom equivalence data.
struct MagneticDataset {
  int uni_number = 0; // 1..1651
  magnetic::MagneticType msg_type = magnetic::MagneticType::type_i;
  // Hall number of the family (types I–III) or maximal (type IV) space group.
  int hall_number = 0;
  SiteTensor tensor_rank = SiteTensor::collinear;

  // Magnetic operations of the input cell.
  MagneticSymmetryOperations operations;
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
// tensors) and standardize it. Port of spglib.c spg_get_magnetic_dataset (3D
// path), wiring spin::operations_with_site_tensors (with time reversal) ->
// magnetic::identify_magnetic_spacegroup_type -> spin::idealized_cell ->
// magnetic::transform_cell. `is_axial` selects axial-vector transformation of
// the rank-1 tensors; `mag_symprec` (std::nullopt -> symprec) is the moment
// tolerance. Errors with e_magnetic_symmetry_search_failed /
// e_cell_standardization_failed on failure.
[[nodiscard]] Result<MagneticDataset>
get_magnetic_dataset(MagneticCell const &cell, bool is_axial,
                     double symprec = kDefaultSymprec,
                     AngleTolerance angle_tolerance = std::nullopt,
                     std::optional<double> mag_symprec = std::nullopt);

} // namespace spglib
