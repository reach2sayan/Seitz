#pragma once

#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <span>
#include <string_view>
#include <utility>

namespace spglib::group {

class SpaceGroup;

// A single Wyckoff position of a space group, as a queryable object (PyXtal's
// Wyckoff_position). Constructed only by SpaceGroup, which owns it. Holds the
// position's multiplicity, Wyckoff letter, site-symmetry symbol, the
// site-symmetry group operations (the stabilizer of a generic point of the
// position), and enough information to expand a free coordinate into the full
// crystallographic orbit.
//
// The data is backed by the compile-time site-symmetry database (decoded via
// data::wyckoff_coordinate / site_symmetry_symbol) and the Hall setting's
// conventional operations (data::operations_from_database). Nothing here is
// lazy — a SpaceGroup builds all of its positions up front.
class WyckoffPosition {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Wyckoff letter index, 0 = 'a'.
  [[nodiscard]] int letter_index() const noexcept { return letter_; }
  // Wyckoff letter as a character (a..z, then A.. for the rare > 26 cases).
  [[nodiscard]] char letter() const noexcept;

  [[nodiscard]] std::string_view site_symmetry() const noexcept {
    return symbol_;
  }

  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  // The site-symmetry group: the operations that fix a generic point of the
  // position. Its order times the multiplicity equals the number of
  // conventional operations of the space group.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return site_symmetry_;
  }

  // The coset representatives that generate the orbit (one per orbit point), and
  // the projector onto the position's canonical coordinate. Together they let a
  // caller expand the orbit with custom folding — the layer-crystal generator
  // uses them to build an orbit without wrapping the non-periodic c axis (which
  // get_all_positions folds, breaking c-flipping operations in a layer cell).
  [[nodiscard]] std::span<SymmetryOperation const>
  orbit_operations() const noexcept {
    return orbit_ops_;
  }
  [[nodiscard]] Vector3d canonical_coordinate(Vector3d const &xyz) const {
    return repr_rotation_.cast<double>() * xyz + repr_translation_;
  }

  // Expand a single free coordinate into the full symmetry orbit (one row per
  // image, folded into the unit cell). `xyz` is first projected onto the
  // position's canonical coordinate, so a point that is only approximately on
  // the position still yields the exact orbit. Row count equals multiplicity()
  // for a generic `xyz` (fewer if `xyz` lands on a higher-symmetry coincidence).
  // PyXtal Wyckoff_position.get_all_positions analogue.
  [[nodiscard]] Positions get_all_positions(Vector3d const &xyz) const;

private:
  friend class SpaceGroup;
  WyckoffPosition(int multiplicity, int letter, int dof,
                  std::string_view symbol, Matrix3i repr_rotation,
                  Vector3d repr_translation, SymmetryOperations orbit_ops,
                  SymmetryOperations site_symmetry)
      : multiplicity_(multiplicity), letter_(letter), dof_(dof),
        symbol_(symbol), repr_rotation_(std::move(repr_rotation)),
        repr_translation_(std::move(repr_translation)),
        orbit_ops_(std::move(orbit_ops)),
        site_symmetry_(std::move(site_symmetry)) {}

  int multiplicity_ = 0;
  int letter_ = 0;
  int dof_ = 0;
  std::string_view symbol_;
  // Representative coordinate operator (projector onto the canonical Wyckoff
  // coordinate): x_canonical = repr_rotation_ . x + repr_translation_.
  Matrix3i repr_rotation_{Matrix3i::Identity()};
  Vector3d repr_translation_{Vector3d::Zero()};
  // Coset representatives generating the orbit (one per orbit point) and the
  // site-symmetry stabilizer, both partitioned from the conventional operations.
  SymmetryOperations orbit_ops_;
  SymmetryOperations site_symmetry_;
};

} // namespace spglib::group
