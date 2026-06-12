#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace spglib::group {

class RodGroup;

// A single Wyckoff position of a rod group (1D-periodic subperiodic group). The
// rod counterpart of group::PointWyckoff: a class of points sharing an orbit
// type, here under a group that carries translations and is periodic along one
// axis. Its locus is an affine subspace `origin + span(basis columns)` of the
// fractional cell; `sample()` turns free parameters into a point on it and
// `orbit_operations()` (coset reps, with translations) expands that point into
// the orbit. Constructed only by RodGroup.
class RodWyckoff {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Free coordinates of the locus (0..3). For the general position this is 3.
  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  [[nodiscard]] char letter() const noexcept { return letter_; }

  // The site-symmetry group (the operations fixing a generic point of the
  // locus). Its order times the multiplicity equals the rod group's order.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return site_symmetry_;
  }

  // Coset representatives generating the orbit (one per orbit point), carrying
  // the operations' translations.
  [[nodiscard]] std::span<SymmetryOperation const>
  orbit_operations() const noexcept {
    return orbit_ops_;
  }

  // A point on the locus (fractional coordinate) from `degrees_of_freedom()`
  // free parameters: origin + sum_i params[i] * basis.col(i).
  [[nodiscard]] Vector3d sample(std::span<double const> params) const;

  // The i-th free direction of the locus (fractional), 0 <= i < dof. The basis
  // is axis-separated: a direction is either the periodic axis (so a free
  // coordinate spanning the repeat) or lies in the aperiodic cross-section — the
  // generator uses this to decide how to sample each free coordinate.
  [[nodiscard]] Vector3d free_direction(int i) const {
    return locus_basis_.col(i);
  }

private:
  friend class RodGroup;
  RodWyckoff(int multiplicity, int dof, char letter, Vector3d locus_origin,
             Matrix3d locus_basis, SymmetryOperations orbit_ops,
             SymmetryOperations site_symmetry)
      : multiplicity_(multiplicity), dof_(dof), letter_(letter),
        locus_origin_(std::move(locus_origin)),
        locus_basis_(std::move(locus_basis)), orbit_ops_(std::move(orbit_ops)),
        site_symmetry_(std::move(site_symmetry)) {}

  int multiplicity_ = 0;
  int dof_ = 0;
  char letter_ = 'a';
  Vector3d locus_origin_{Vector3d::Zero()};
  // Columns 0..dof_-1 span the locus directions; the rest are unused.
  Matrix3d locus_basis_{Matrix3d::Zero()};
  SymmetryOperations orbit_ops_;
  SymmetryOperations site_symmetry_;
};

// A rod group (1..75) as a standalone, structure-free object: it owns its
// symmetry operations (decoded from the generated rod tables) and its Wyckoff
// positions, derived in-house from those operations. The 1D-periodic sibling of
// group::SpaceGroup (3D), SpaceGroup::from_layer_* (2D) and group::PointGroup
// (0D).
//
// The Wyckoff positions (general AND special) are derived from the operations by
// the affine fixed-locus arrangement (rod_group.cpp) — the translation-bearing
// generalisation of group::PointGroup's linear subspace arrangement, with the
// periodic axis handled modulo the rod lattice: each operation's fixed locus is
// the affine solution of (R - I) p = -t + n*e_axis, the loci are closed under
// intersection, and one Wyckoff position is emitted per orbit of loci.
class RodGroup {
public:
  // Build by rod-group number (1..75). Requires the generated rod tables
  // (data/rod_group_tables.hpp via tools/transcribe_rod_groups.py).
  [[nodiscard]] static Result<RodGroup> from_number(int number);

  [[nodiscard]] int number() const noexcept { return number_; }
  [[nodiscard]] std::string_view symbol() const noexcept { return symbol_; }

  // The single periodic axis (0=a, 1=b, 2=c); orbit expansion folds only this
  // component (the other two are aperiodic / vacuum directions).
  [[nodiscard]] int periodic_axis() const noexcept { return periodic_axis_; }

  [[nodiscard]] int order() const noexcept {
    return static_cast<int>(operations_.size());
  }

  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return operations_;
  }

  // The Wyckoff positions, ordered by ascending multiplicity ('a' = the most
  // special); the last is the general position.
  [[nodiscard]] std::span<RodWyckoff const> wyckoffs() const noexcept {
    return positions_;
  }

private:
  RodGroup() = default;

  int number_ = 0;
  std::string_view symbol_;
  int periodic_axis_ = 2;
  SymmetryOperations operations_;
  std::vector<RodWyckoff> positions_;
};

} // namespace spglib::group
