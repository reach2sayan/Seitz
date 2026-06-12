#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace spglib::group {

class PointGroup;

// A single Wyckoff position of a crystallographic point group: a class of
// points sharing an orbit type under the point group acting about the origin.
// Unlike the space-group WyckoffPosition there is no translation part and no
// site-symmetry database — everything is derived from the point group's
// operations (the arrangement of fixed subspaces: the origin, the rotation
// axes, the mirror planes, and their intersections).
//
// The position is parameterised by a `degrees_of_freedom()`-dimensional locus
// (its fixed subspace): `sample()` turns free parameters into a point on it, and
// `orbit_operations()` expands that point into the full orbit. Constructed only
// by PointGroup, which owns it.
class PointWyckoff {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Free coordinates of the locus: 0 (the origin), 1 (a point on an axis),
  // 2 (a point on a plane), or 3 (the general position).
  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  // Wyckoff letter, 'a' = the most special (smallest multiplicity).
  [[nodiscard]] char letter() const noexcept { return letter_; }

  // The site-symmetry group: the operations that fix a generic point of the
  // position. Its order times the multiplicity equals the point group's order.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return site_symmetry_;
  }

  // The coset representatives that generate the orbit (one per orbit point).
  // Applying each to a point from sample() yields the full orbit.
  [[nodiscard]] std::span<SymmetryOperation const>
  orbit_operations() const noexcept {
    return orbit_ops_;
  }

  // A point on the position (fractional coordinate in the conventional basis)
  // from `degrees_of_freedom()` free parameters. Extra parameters are ignored;
  // a 0-DOF position returns the fixed point (origin) regardless of `params`.
  [[nodiscard]] Vector3d sample(std::span<double const> params) const;

private:
  friend class PointGroup;
  PointWyckoff(int multiplicity, int dof, char letter, Matrix3d locus_basis,
               SymmetryOperations orbit_ops, SymmetryOperations site_symmetry)
      : multiplicity_(multiplicity), dof_(dof), letter_(letter),
        locus_basis_(std::move(locus_basis)), orbit_ops_(std::move(orbit_ops)),
        site_symmetry_(std::move(site_symmetry)) {}

  int multiplicity_ = 0;
  int dof_ = 0;
  char letter_ = 'a';
  // Columns 0..dof_-1 are an orthonormal basis of the locus (fractional basis);
  // the remaining columns are unused.
  Matrix3d locus_basis_{Matrix3d::Zero()};
  SymmetryOperations orbit_ops_;
  SymmetryOperations site_symmetry_;
};

// A crystallographic point group (1..32) as a standalone, structure-free
// object: it owns its symmetry operations (pure rotations / rotoinversions
// about the origin, expressed as integer matrices in the conventional basis)
// and its Wyckoff positions as first-class objects.
//
// The operations are taken from a representative symmorphic space group of the
// point group (the distinct rotation parts of its conventional operations); the
// Wyckoff positions are then derived from those operations. This is the 0D
// counterpart of SpaceGroup (3D) and SpaceGroup::from_layer_* (2D) for
// crystal/cluster generation.
class PointGroup {
public:
  // Build by international point-group number (1..32, the spglib pointgroup
  // numbering: 1, -1, 2, m, 2/m, 222, ...).
  [[nodiscard]] static Result<PointGroup> from_number(int number);

  [[nodiscard]] int number() const noexcept { return number_; }
  [[nodiscard]] std::string_view symbol() const noexcept { return symbol_; }
  [[nodiscard]] std::string_view schoenflies() const noexcept {
    return schoenflies_;
  }

  // A space group whose point group is this one (a symmorphic P group); its
  // crystal system fixes the metric used when realising a cluster in Cartesian
  // space (generate::random_cluster).
  [[nodiscard]] int representative_spacegroup() const noexcept {
    return representative_spacegroup_;
  }

  // The order of the point group.
  [[nodiscard]] int order() const noexcept {
    return static_cast<int>(operations_.size());
  }

  // The symmetry operations (integer rotations, zero translation).
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return operations_;
  }

  // The Wyckoff positions, ordered by ascending multiplicity ('a' = the most
  // special); the last is the general position.
  [[nodiscard]] std::span<PointWyckoff const> wyckoffs() const noexcept {
    return positions_;
  }

private:
  PointGroup() = default;

  int number_ = 0;
  std::string_view symbol_;
  std::string_view schoenflies_;
  int representative_spacegroup_ = 0;
  SymmetryOperations operations_;
  std::vector<PointWyckoff> positions_;
};

} // namespace spglib::group
