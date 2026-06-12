#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>

#include <optional>
#include <vector>

// Magnetic symmetry from per-site tensors (3D path). Given the non-magnetic
// (spatial) symmetry of a cell and its per-site magnetic moments, this
// filters/splits those operations into the magnetic symmetry operations that
// also preserve the moments (optionally up to time reversal).
namespace spglib::spin {

// The result of the magnetic symmetry search, in the *input* cell's basis.
struct MagneticSymmetrySearch {
  MagneticSymmetryOperations operations;
  // equivalent_atoms[i] = representative atom of i's magnetic orbit.
  std::vector<int> equivalent_atoms;
  // permutations[p * size + i] = image of atom i under operation p.
  std::vector<int> permutations;
  // Primitive lattice implied by the magnetic pure translations.
  Matrix3d primitive_lattice{Matrix3d::Identity()};
};

// Magnetic symmetry operations of `mcell` consistent with its site tensors.
//   - with_time_reversal: when true, anti-operations (moment-reversing) are
//     allowed and the result is the magnetic family space group; when false,
//     only ordinary operations are kept (the maximal space subgroup).
//   - is_axial: site tensors of rank 1 transform as axial vectors
//     (v' = |det R| R v) rather than ordinary vectors (v' = R v). Collinear
//     scalars pick up a |det R| factor when is_axial.
//   - mag_symprec: tolerance for comparing moments; std::nullopt uses `symprec`.
// `sym_nonspin` is the spatial symmetry (e.g. from symmetry::find_symmetry).
// Errors with e_magnetic_symmetry_search_failed when the orbits or the
// primitive lattice cannot be resolved. For the object-oriented entry point that
// derives the spatial symmetry for you and memoizes, prefer
// spglib::analysis::MagneticSymmetryAnalyzer::symmetry_search().
[[nodiscard]] Result<MagneticSymmetrySearch>
operations_with_site_tensors(SymmetryOperations const &sym_nonspin,
                             MagneticCell const &mcell, bool with_time_reversal,
                             bool is_axial, double symprec,
                             AngleTolerance angle_tolerance = std::nullopt,
                             std::optional<double> mag_symprec = std::nullopt);

// The pure translations of a magnetic symmetry: the translations whose rotation
// is the identity and which are not time-reversal operations (the type-IV
// anti-translations are excluded). Includes the zero translation. For the
// object-oriented form, prefer
// spglib::analysis::MagneticOperationSet::pure_translations().
[[nodiscard]] std::vector<Vector3d>
collect_pure_translations(MagneticSymmetryOperations const &operations);

// The idealized magnetic cell: each atom's position and site tensor is replaced
// by the average, over all magnetic operations, of the operation applied to the
// atom that maps onto it. `permutations[p * size + i]` is the image of atom i
// under operation p (as returned by operations_with_site_tensors).
[[nodiscard]] MagneticCell
idealized_cell(std::vector<int> const &permutations, MagneticCell const &mcell,
               MagneticSymmetryOperations const &operations,
               bool with_time_reversal, bool is_axial);

} // namespace spglib::spin
