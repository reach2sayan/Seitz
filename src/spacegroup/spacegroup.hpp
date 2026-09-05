#pragma once

#include "core/testable.hpp"
#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/spacegroup_match.hpp>

#include "symmetry/primitive.hpp"

#include <optional>

namespace seitz::spacegroup {

// Matches a primitive cell against the Hall-symbol database. The family is a
// compile-time parameter: the layer path draws its candidates from the layer
// catalog, tries only the two c-preserving orthorhombic axis choices, and skips
// the 3D-only representative-Hall refinement -- all `if constexpr` branches.
//
// Non-owning: `primitive` must outlive the matcher.
template <GroupFamily F> class SEITZ_TESTABLE SpacegroupMatcher {
public:
  // An unset `setting` searches every Hall setting of the family; a set one
  // restricts the search to it. The tolerance is the one the primitive cell was
  // found at.
  SpacegroupMatcher(symmetry::Primitive const &primitive,
                    std::optional<HallNumber> setting) noexcept
      : primitive_{primitive}, setting_{setting} {}

  // The matched setting. Errors with e_spacegroup_search_failed.
  [[nodiscard]] Result<SpacegroupMatch> search() const;

  // Try one Hall setting against operations already in the bravais /
  // conventional setting. Returns the origin shift that reproduces that
  // setting's database operations, std::nullopt otherwise.
  [[nodiscard]] static std::optional<Vector3d>
  match_hall(Matrix3d const &bravais_lattice, HallNumber hall,
             data::Centering centering, Operations const &symmetry,
             double symprec);

private:
  symmetry::Primitive const &primitive_;
  std::optional<HallNumber> setting_;
};

// search() lives in spacegroup.cpp and match_hall() in hall_symbol.cpp, so the
// members are instantiated one at a time rather than by whole-class
// instantiation.
extern template Result<SpacegroupMatch>
SpacegroupMatcher<GroupFamily::space>::search() const;
extern template Result<SpacegroupMatch>
SpacegroupMatcher<GroupFamily::layer>::search() const;
extern template std::optional<Vector3d>
SpacegroupMatcher<GroupFamily::space>::match_hall(Matrix3d const &, HallNumber,
                                                  data::Centering,
                                                  Operations const &, double);
extern template std::optional<Vector3d>
SpacegroupMatcher<GroupFamily::layer>::match_hall(Matrix3d const &, HallNumber,
                                                  data::Centering,
                                                  Operations const &, double);

// Determine the space group of a primitive cell given by a lattice and its
// symmetry operations (a single notional atom at the origin). The magnetic
// path needs this without the Niggli reduction that
// OperationSet::spacegroup applies.
[[nodiscard]] Result<SpacegroupMatch>
search_spacegroup_with_symmetry(Operations const &operations,
                                Matrix3d const &prim_lattice, double symprec);

} // namespace seitz::spacegroup
