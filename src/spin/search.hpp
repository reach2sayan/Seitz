#pragma once

#include "core/testable.hpp"
#include <seitz/core/error.hpp>
#include <seitz/core/magnetic_cell.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/mdspan.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

// Magnetic symmetry from per-site tensors (3D path). Given the non-magnetic
// (spatial) symmetry of a cell and its per-site magnetic moments, this
// filters/splits those operations into the magnetic symmetry operations that
// also preserve the moments (optionally up to time reversal).
namespace seitz::spin {

// The result of the magnetic symmetry search, in the *input* cell's basis.
class MagneticSymmetrySearch {
public:
  MagneticOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;
  // Primitive lattice implied by the magnetic pure translations.
  Matrix3d primitive_lattice{Matrix3d::Identity()};

  MagneticSymmetrySearch(MagneticOperations ops, std::vector<int> equivalent,
                         std::vector<int> permutations, Matrix3d prim_lattice)
      : operations{std::move(ops)}, equivalent_atoms{std::move(equivalent)},
        primitive_lattice{std::move(prim_lattice)},
        permutations_{std::move(permutations)} {}

  // permutations()[p, i] = image of atom i under operation p: a
  // (#operations x #atoms) view over the private flat buffer.
  [[nodiscard]] md::matrix_view<int const> permutations() const noexcept {
    return md::matrix_view<int const>(
        permutations_.data(), static_cast<Index>(operations.size()),
        static_cast<Index>(equivalent_atoms.size()));
  }

private:
  // One row of equivalent_atoms.size() images per operation.
  std::vector<int> permutations_;
};

// The magnetic symmetry search over one cell: which of the cell's spatial
// operations survive once the site tensors must be preserved too, and the
// idealized cell that follows. The tensor kind (axial vs polar) is carried by
// the cell; time reversal is a compile-time parameter of the search, since it
// selects between the magnetic family group and its maximal space subgroup.
//
// Non-owning: `cell` and `spatial` must outlive the search.
class SEITZ_TESTABLE SpinSearch {
public:
  SpinSearch(MagneticCell const &cell, Operations const &spatial,
             MagneticTolerance const &tol) noexcept
      : cell_{cell}, spatial_{spatial}, tol_{tol} {}

  // The magnetic operations, orbits, permutations and primitive lattice.
  // `TimeReversal::on` allows anti-operations (moment-reversing) and yields the
  // magnetic family space group; `off` keeps only ordinary operations. Errors
  // with e_magnetic_symmetry_search_failed when the orbits or the primitive
  // lattice cannot be resolved.
  template <TimeReversal TR>
  [[nodiscard]] Result<MagneticSymmetrySearch> operations() const;

  // The idealized magnetic cell: each atom's position and site tensor replaced
  // by the average, over all magnetic operations, of the operation applied to
  // the atom that maps onto it.
  template <TimeReversal TR>
  [[nodiscard]] MagneticCell
  idealized(MagneticSymmetrySearch const &search) const;

private:
  MagneticCell const &cell_;
  Operations const &spatial_;
  MagneticTolerance tol_;
};

extern template Result<MagneticSymmetrySearch>
SpinSearch::operations<TimeReversal::on>() const;
extern template Result<MagneticSymmetrySearch>
SpinSearch::operations<TimeReversal::off>() const;
extern template MagneticCell
SpinSearch::idealized<TimeReversal::on>(MagneticSymmetrySearch const &) const;
extern template MagneticCell
SpinSearch::idealized<TimeReversal::off>(MagneticSymmetrySearch const &) const;

} // namespace seitz::spin
