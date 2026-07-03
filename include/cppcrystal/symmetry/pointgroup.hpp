#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <span>

namespace cppcrystal::symmetry {

// Point-group data (symbol, Schoenflies, holohedry, Laue) for number 1..32;
// number 0 / out of range gives an empty PointGroup.
[[nodiscard]] PointGroup pointgroup_by_number(int number) noexcept;

// Identify the point-group number (1..32) from a set of rotations (duplicates
// allowed); 0 if none matches.
[[nodiscard]] int
identify_pointgroup_number(std::span<Matrix3i const> rotations) noexcept;

// Point group together with the integer change-of-basis matrix that brings the
// rotations into the conventional setting (columns are the chosen axes). With
// `aperiodic_axis` set (a layer cell), the cubic point groups are rejected and
// the conventional axes are chosen so the aperiodic axis becomes c.
struct PointgroupTransform {
  PointGroup pointgroup;
  Matrix3i transformation{Matrix3i::Zero()};
};

[[nodiscard]] Result<PointgroupTransform>
get_pointgroup(std::span<Matrix3i const> rotations,
               std::optional<int> aperiodic_axis = std::nullopt);

// Convenience: extract the (de-duplicated by value is done internally) rotation
// parts of a set of symmetry operations. For the object-oriented form, prefer
// cppcrystal::analysis::OperationSet::rotations().
[[nodiscard]] std::vector<Matrix3i> rotations_of(SymmetryOperations const &ops);

// Crystallographic type of a single rotation, in a lattice basis. The sign of
// the corresponding order (below) is proper (+) vs improper (−); the magnitude
// is the rotation order.
enum class RotationType {
  rotoinversion_6, // -6
  rotoinversion_4, // -4
  rotoinversion_3, // -3
  mirror,          // -2
  inversion,       // -1
  identity,        //  1
  rotation_2,      //  2
  rotation_3,      //  3
  rotation_4,      //  4
  rotation_6,      //  6
};

// Classify a single integer rotation by determinant and trace; std::nullopt if
// it is not a crystallographic rotation.
[[nodiscard]] std::optional<RotationType> rotation_type(Matrix3i const &rot) noexcept;

// Signed crystallographic order of a rotation (International Tables convention):
// +n for an n-fold proper rotation (n ∈ {1,2,3,4,6}); −1 inversion, −2 mirror,
// −3/−4/−6 rotoinversion. 0 if the matrix is not a crystallographic rotation.
[[nodiscard]] int rotation_order(Matrix3i const &rot) noexcept;

// The invariant axis of a rotation's proper part, expressed in the lattice
// basis: the rotation axis, or the plane normal for a mirror. std::nullopt for
// the identity and inversion (no unique axis) or a non-crystallographic matrix.
[[nodiscard]] std::optional<Vector3i> rotation_axis(Matrix3i const &rot);

} // namespace cppcrystal::symmetry
