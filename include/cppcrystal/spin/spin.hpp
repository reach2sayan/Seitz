#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/mdspan.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

// Magnetic symmetry from per-site tensors (3D path). Given the non-magnetic
// (spatial) symmetry of a cell and its per-site magnetic moments, this
// filters/splits those operations into the magnetic symmetry operations that
// also preserve the moments (optionally up to time reversal).
namespace cppcrystal::spin {

// The result of the magnetic symmetry search, in the *input* cell's basis.
class MagneticSymmetrySearch {
public:
  MagneticSymmetryOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;
  // Primitive lattice implied by the magnetic pure translations.
  Matrix3d primitive_lattice{Matrix3d::Identity()};

  MagneticSymmetrySearch(MagneticSymmetryOperations ops,
                         std::vector<int> equivalent,
                         std::vector<int> permutations, Matrix3d prim_lattice)
      : operations(std::move(ops)), equivalent_atoms(std::move(equivalent)),
        primitive_lattice(std::move(prim_lattice)),
        permutations_(std::move(permutations)) {}

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

// Magnetic symmetry operations of `mcell` consistent with its site tensors.
//   - time_reversal: `on` allows anti-operations (moment-reversing) and yields
//     the magnetic family space group; `off` keeps only ordinary operations
//     (the maximal space subgroup).
//   - the tensor kind (axial vs polar) is carried by `mcell`.
//   - tol.moment: tolerance for comparing moments; unset uses tol.symprec.
// `sym_nonspin` is the spatial symmetry (e.g. from symmetry::find_symmetry).
// Errors with e_magnetic_symmetry_search_failed when the orbits or the
// primitive lattice cannot be resolved. For the object-oriented entry point that
// derives the spatial symmetry for you and memoizes, prefer
// cppcrystal::analysis::MagneticSymmetryAnalyzer::symmetry_search().
[[nodiscard]] Result<MagneticSymmetrySearch>
operations_with_site_tensors(SymmetryOperations const &sym_nonspin,
                             MagneticCell const &mcell,
                             TimeReversal time_reversal,
                             MagneticTolerance const &tol);

// The pure translations of a magnetic symmetry: the translations whose rotation
// is the identity and which are not time-reversal operations (the type-IV
// anti-translations are excluded). Includes the zero translation. For the
// object-oriented form, prefer
// cppcrystal::analysis::MagneticOperationSet::pure_translations().
[[nodiscard]] std::vector<Vector3d>
collect_pure_translations(MagneticSymmetryOperations const &operations);

// The idealized magnetic cell: each atom's position and site tensor is replaced
// by the average, over all magnetic operations, of the operation applied to the
// atom that maps onto it (per `search`'s operations and permutations).
[[nodiscard]] MagneticCell idealized_cell(MagneticSymmetrySearch const &search,
                                          MagneticCell const &mcell,
                                          TimeReversal time_reversal);

} // namespace cppcrystal::spin
