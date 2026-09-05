#pragma once

#include "core/testable.hpp"
#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/point_group.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#include <optional>
#include <span>

namespace seitz::symmetry {

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
[[nodiscard]] SEITZ_TESTABLE Result<PointgroupTransform>
identify_point_group(std::span<Matrix3i const> rotations,
                     std::optional<int> layer_axis = std::nullopt);

extern template Result<PointgroupTransform>
    identify_point_group<GroupFamily::space>(std::span<Matrix3i const>,
                                             std::optional<int>);
extern template Result<PointgroupTransform>
    identify_point_group<GroupFamily::layer>(std::span<Matrix3i const>,
                                             std::optional<int>);

} // namespace seitz::symmetry
