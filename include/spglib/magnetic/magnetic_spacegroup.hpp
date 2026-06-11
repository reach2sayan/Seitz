#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/types.hpp>

// Magnetic space-group determination. Port of magnetic_spacegroup.c (3D path):
// given a lattice and a set of magnetic symmetry operations, identify the
// magnetic space group (UNI number) and the transformation to its standardized
// setting. Builds on spacegroup::{primitive_symmetry, search_spacegroup_with_
// symmetry}, refine::{find_similar_bravais_lattice, conventional_lattice}, the
// msg_database, and spin::collect_pure_translations.
namespace spglib::magnetic {

// Construction type of a magnetic space group (Bärnighausen / BNS types I–IV).
// Replaces spglib's bare int `type` 1..4. Numeric values match the database
// `type` field so the two can be bridged at the boundary.
enum class MagneticType { type_i = 1, type_ii = 2, type_iii = 3, type_iv = 4 };

// The identification result — port of magnetic_spacegroup.h MagneticTypeIdentification
// (the non-cell part). `transformation_matrix` and `origin_shift` map the input
// setting to the standardized one; `std_rotation_matrix` is the rigid rotation
// to the idealized standardized lattice.
struct MagneticTypeIdentification {
  int uni_number = 0; // 1..1651
  MagneticType msg_type = MagneticType::type_i;
  int hall_number = 0; // family (types I–III) or maximal (type IV) space group

  Matrix3d transformation_matrix{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};
  Matrix3d std_rotation_matrix{Matrix3d::Identity()};
};

// Identify the magnetic space-group type of `magnetic_symmetry` in `lattice`.
// Port of msg_identify_magnetic_space_group_type. Errors with
// e_magnetic_symmetry_search_failed when no UNI number matches.
[[nodiscard]] Result<MagneticTypeIdentification>
identify_magnetic_spacegroup_type(
    Matrix3d const &lattice,
    MagneticSymmetryOperations const &magnetic_symmetry, double symprec);

// Transform an (idealized) magnetic cell to the standardized setting given the
// transformation (transformation_matrix, origin_shift) and the rigid rotation
// to the idealized lattice. The cell size may change (primitive -> conventional
// centering). Site tensors are Cartesian, so only the rigid rotation acts on
// them (rank-1) or nothing (rank-0). Port of msg_get_transformed_cell. Errors
// with e_cell_standardization_failed.
[[nodiscard]] Result<MagneticCell>
transform_cell(MagneticCell const &mcell, Matrix3d const &transformation_matrix,
               Vector3d const &origin_shift, Matrix3d const &rigid_rotation,
               MagneticSymmetryOperations const &magnetic_symmetry,
               double symprec);

} // namespace spglib::magnetic
