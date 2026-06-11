#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/symmetry/primitive.hpp>

namespace spglib::spacegroup {

// The determined space group of a (primitive) cell: the database metadata
// (number, Hall number, symbols, point group, centering) together with the
// conventional/bravais lattice it was matched in and the origin shift that
// aligns the cell's operations with the Hall-symbol database. Value type port
// of spacegroup.c's Spacegroup struct (the string fields live in `type`, which
// is a view into the static database tables).
struct Spacegroup {
  data::SpacegroupType type;
  Matrix3d bravais_lattice{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};

  [[nodiscard]] int number() const noexcept { return type.number; }
  [[nodiscard]] int hall_number() const noexcept { return type.hall_number; }
};

// Determine the space group of a primitive cell. Port of spacegroup.c
// spa_search_spacegroup (3D space-group path). `hall_number == 0` searches all
// 230 space groups; a non-zero value restricts the search to that Hall setting.
//
// Errors with e_spacegroup_search_failed when no Hall setting matches.
[[nodiscard]] Result<Spacegroup>
search_spacegroup(symmetry::Primitive const &primitive, int hall_number,
                  double symprec, AngleTolerance angle_tolerance);

// Whether the `lattice` given to spacegroup_type_from_symmetry is the
// conventional cell (recover the primitive setting via the t_mat implied by the
// operations) or already a primitive cell. Replaces spglib's
// `transform_lattice_by_tmat` int flag — it is a compile-time template
// argument, so the branch is resolved at compile time.
enum class LatticeSetting { conventional, primitive };

// Determine the space group directly from a set of symmetry operations (no
// atomic positions). Port of spglib.c get_hall_number_from_symmetry: build the
// primitive symmetry (symmetry::primitive_symmetry), Niggli-reduce the implied
// primitive lattice, and match it against the Hall database.
// `LatticeSetting::conventional` mirrors spg_get_spacegroup_type_from_symmetry;
// `LatticeSetting::primitive` (with an identity lattice) mirrors
// spg_get_hall_number_from_symmetry. Errors with e_spacegroup_search_failed.
template <LatticeSetting Setting = LatticeSetting::conventional>
[[nodiscard]] Result<Spacegroup>
spacegroup_type_from_symmetry(SymmetryOperations const &operations,
                              Matrix3d const &lattice, double symprec);

// Determine the space group of a primitive cell given by a lattice and its
// symmetry operations (a single notional atom at the origin). Port of
// spacegroup.c spa_search_spacegroup_with_symmetry.
[[nodiscard]] Result<Spacegroup>
search_spacegroup_with_symmetry(SymmetryOperations const &operations,
                                Matrix3d const &prim_lattice, double symprec);

} // namespace spglib::spacegroup
