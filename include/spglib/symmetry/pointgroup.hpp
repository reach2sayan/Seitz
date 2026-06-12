#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/point_group.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <optional>
#include <span>

namespace spglib::symmetry {

// Point-group data (symbol, Schoenflies, holohedry, Laue) for number 1..32;
// number 0 / out of range gives an empty PointGroup. Port of
// ptg_get_pointgroup.
[[nodiscard]] PointGroup pointgroup_by_number(int number) noexcept;

// Identify the point-group number (1..32) from a set of rotations (duplicates
// allowed); 0 if none matches. Port of get_pointgroup_number_by_rotations.
[[nodiscard]] int
identify_pointgroup_number(std::span<Matrix3i const> rotations) noexcept;

// Point group together with the integer change-of-basis matrix that brings the
// rotations into the conventional setting (columns are the chosen axes). Port
// of ptg_get_transformation_matrix. With `aperiodic_axis` set (a layer cell),
// the cubic point groups are rejected and the conventional axes are chosen so
// the aperiodic axis becomes c.
struct PointgroupTransform {
  PointGroup pointgroup;
  Matrix3i transformation{Matrix3i::Zero()};
};

[[nodiscard]] Result<PointgroupTransform>
get_pointgroup(std::span<Matrix3i const> rotations,
               std::optional<int> aperiodic_axis = std::nullopt);

// Convenience: extract the (de-duplicated by value is done internally) rotation
// parts of a set of symmetry operations. For the object-oriented form, prefer
// spglib::analysis::OperationSet::rotations().
[[nodiscard]] std::vector<Matrix3i> rotations_of(SymmetryOperations const &ops);

} // namespace spglib::symmetry
