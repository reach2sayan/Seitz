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

// Convenience: the rotation parts of a set of symmetry operations, one per
// operation in order (de-duplication happens in get_pointgroup). For the object-oriented form, prefer
// cppcrystal::analysis::OperationSet::rotations().
[[nodiscard]] std::vector<Matrix3i> rotations_of(SymmetryOperations const &ops);

} // namespace cppcrystal::symmetry
