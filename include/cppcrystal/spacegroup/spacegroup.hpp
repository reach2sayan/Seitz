#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/symmetry/primitive.hpp>

// For the object-oriented entry points wrapping a bare operation set, prefer
// cppcrystal::analysis::OperationSet::{spacegroup_type, search_spacegroup}.
namespace cppcrystal::spacegroup {

// The determined space group of a (primitive) cell: the database metadata
// (number, Hall number, symbols, point group, centering) together with the
// conventional/bravais lattice it was matched in and the origin shift that
// aligns the cell's operations with the Hall-symbol database. The string fields
// live in `type`, which is a view into the static database tables.
struct Spacegroup {
  data::SpacegroupType type;
  Matrix3d bravais_lattice{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};

  [[nodiscard]] int number() const noexcept { return type.number; }
  [[nodiscard]] int hall_number() const noexcept { return type.hall_number; }
};

// Determine the space group of a primitive cell. `hall_number == 0` searches all
// 230 space groups; a non-zero value restricts the search to that Hall setting.
[[nodiscard]] Result<Spacegroup>
search_spacegroup(symmetry::Primitive const &primitive, int hall_number,
                  Tolerance const &tol);

// Whether the `lattice` given to spacegroup_type_from_symmetry is the
// conventional cell (recover the primitive setting via the t_mat implied by the
// operations) or already a primitive cell. A compile-time template argument, so
// the branch is resolved at compile time.
enum class LatticeSetting { conventional, primitive };

// Determine the space group directly from a set of symmetry operations (no
// atomic positions): build the primitive symmetry, Niggli-reduce the implied
// primitive lattice, and match it against the Hall database.
// `LatticeSetting::primitive` is used with an identity lattice. Errors with
// e_spacegroup_search_failed.
template <LatticeSetting Setting = LatticeSetting::conventional>
[[nodiscard]] Result<Spacegroup>
spacegroup_type_from_symmetry(SymmetryOperations const &operations,
                              Matrix3d const &lattice, double symprec);

// Determine the space group of a primitive cell given by a lattice and its
// symmetry operations (a single notional atom at the origin).
[[nodiscard]] Result<Spacegroup>
search_spacegroup_with_symmetry(SymmetryOperations const &operations,
                                Matrix3d const &prim_lattice, double symprec);

} // namespace cppcrystal::spacegroup
