#pragma once

#include "core/testable.hpp"
#include <cppcrystal/analysis/dataset.hpp> // magnetic::MagneticType
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

// Magnetic space-group determination (3D path): given a lattice and a set of
// magnetic symmetry operations, identify the magnetic space group (UNI number)
// and the transformation to its standardized setting.
namespace cppcrystal::magnetic {

// The identification result. `transformation_matrix` and `origin_shift` map the
// input setting to the standardized one; `std_rotation_matrix` is the rigid
// rotation to the idealized standardized lattice.
struct MagneticTypeIdentification {
  UniNumber uni;
  MagneticType msg_type = MagneticType::type_i;
  // Family (types I-III) or maximal (type IV) space group.
  HallNumber hall;

  Matrix3d transformation_matrix{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};
  Matrix3d std_rotation_matrix{Matrix3d::Identity()};
};

// Identifies a magnetic space group from its operations in a given lattice,
// and transforms a cell into the standardized setting that follows. 3D path
// only: magnetic layer groups are not in the database.
//
// Non-owning: `lattice` and `operations` must outlive the identification.
class CPPCRYSTAL_TESTABLE MagneticIdentification {
public:
  MagneticIdentification(Lattice const &lattice,
                         MagneticOperations const &operations,
                         Tolerance const &tol) noexcept
      : lattice_(lattice), operations_(operations), tol_(tol) {}

  // The UNI number, MSG type, and the transformation into the standardized
  // setting. Errors with e_magnetic_symmetry_search_failed when no UNI number
  // matches.
  [[nodiscard]] Result<MagneticTypeIdentification> identify() const;

  // Transform an (idealized) magnetic cell into the standardized setting named
  // by `identification`. The cell size may change (primitive -> conventional
  // centering). Site tensors are Cartesian, so only the rigid rotation acts on
  // them (rank-1) or nothing (rank-0). Errors with
  // e_cell_standardization_failed.
  [[nodiscard]] Result<MagneticCell>
  transform(MagneticCell const &mcell,
            MagneticTypeIdentification const &identification) const;

private:
  Lattice const &lattice_;
  MagneticOperations const &operations_;
  Tolerance tol_;
};

} // namespace cppcrystal::magnetic
