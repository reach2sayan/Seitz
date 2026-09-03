#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/spg_database.hpp>

#include "symmetry/primitive.hpp"

#include <optional>

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

// Matches a primitive cell against the Hall-symbol database. The family is a
// compile-time parameter: the layer path draws its candidates from the layer
// catalog, tries only the two c-preserving orthorhombic axis choices, and skips
// the 3D-only representative-Hall refinement -- all `if constexpr` branches.
//
// Non-owning: `primitive` must outlive the matcher.
template <GroupFamily F> class SpacegroupMatcher {
public:
  // `hall_number == 0` searches every setting of the family; a non-zero value
  // (signed by the database convention) restricts the search to it.
  SpacegroupMatcher(symmetry::Primitive const &primitive,
                    int hall_number) noexcept
      : primitive_(primitive), hall_number_(hall_number) {}

  // The matched setting. Errors with e_spacegroup_search_failed.
  [[nodiscard]] Result<Spacegroup> search(Tolerance const &tol) const;

  // Try one Hall setting against operations already in the bravais /
  // conventional setting. Returns the origin shift that reproduces that
  // setting's database operations, std::nullopt otherwise.
  [[nodiscard]] static std::optional<Vector3d>
  match_hall(Matrix3d const &bravais_lattice, int hall_number,
             data::Centering centering, Operations const &symmetry,
             double symprec);

private:
  symmetry::Primitive const &primitive_;
  int hall_number_;
};

extern template class SpacegroupMatcher<GroupFamily::space>;
extern template class SpacegroupMatcher<GroupFamily::layer>;

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
spacegroup_type_from_symmetry(Operations const &operations,
                              Matrix3d const &lattice, double symprec);

// Determine the space group of a primitive cell given by a lattice and its
// symmetry operations (a single notional atom at the origin).
[[nodiscard]] Result<Spacegroup>
search_spacegroup_with_symmetry(Operations const &operations,
                                Matrix3d const &prim_lattice, double symprec);

} // namespace cppcrystal::spacegroup
