#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/spacegroup_match.hpp>

// Magnetic space-group determination (3D path): given a lattice and a set of
// magnetic symmetry operations, identify the magnetic space group (UNI number)
// and the transformation to its standardized setting.
namespace seitz::magnetic {

// Identifies a magnetic space group from its operations in a given lattice,
// and transforms a cell into the standardized setting that follows. 3D path
// only: magnetic layer groups are not in the database.
//
// Non-owning: `lattice` and `operations` must outlive the identification.
class MagneticIdentification {
public:
  MagneticIdentification(Lattice const &lattice,
                         MagneticOperations const &operations,
                         Tolerance const &tol) noexcept
      : lattice_(lattice), operations_(operations), tol_(tol) {}

  // The UNI number, MSG type, and the transformation into the standardized
  // setting. Errors with e_magnetic_symmetry_search_failed when no UNI number
  // matches.
  [[nodiscard]] Result<MagneticMatch> identify() const;

  // Transform an (idealized) magnetic cell into the standardized setting named
  // by `identification`. The cell size may change (primitive -> conventional
  // centering). Site tensors are Cartesian, so only the rigid rotation acts on
  // them (rank-1) or nothing (rank-0). Errors with
  // e_cell_standardization_failed.
  [[nodiscard]] Result<MagneticCell>
  transform(MagneticCell const &mcell, MagneticMatch const &match) const;

private:
  Lattice const &lattice_;
  MagneticOperations const &operations_;
  Tolerance tol_;
};

} // namespace seitz::magnetic
