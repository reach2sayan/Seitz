#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <span>

namespace cppcrystal::symmetry {


// Point group together with the integer change-of-basis matrix that brings the
// rotations into the conventional setting (columns are the chosen axes). With
// `aperiodic_axis` set (a layer cell), the cubic point groups are rejected and
// the conventional axes are chosen so the aperiodic axis becomes c.
struct PointgroupTransform {
  PointGroup pointgroup;
  Matrix3i transformation{Matrix3i::Zero()};
};

// The family is a compile-time parameter: only the layer path rejects the
// cubic point groups and sorts the axes so the aperiodic one becomes c.
// `layer_axis` is still data — a layer cell's aperiodic axis is whichever of
// the three the input basis put it on — and is ignored for GroupFamily::space.
template <GroupFamily F>
[[nodiscard]] Result<PointgroupTransform>
identify_point_group(std::span<Matrix3i const> rotations,
                     std::optional<int> layer_axis = std::nullopt);

extern template Result<PointgroupTransform>
identify_point_group<GroupFamily::space>(std::span<Matrix3i const>,
                                         std::optional<int>);
extern template Result<PointgroupTransform>
identify_point_group<GroupFamily::layer>(std::span<Matrix3i const>,
                                         std::optional<int>);

} // namespace cppcrystal::symmetry
