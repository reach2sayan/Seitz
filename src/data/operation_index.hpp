#pragma once

#include <cppcrystal/data/spg_database.hpp>

#include "core/matrix_order.hpp"

// The rotation index over a database setting's operations. Private: only the
// Hall-symbol matcher needs it, and RotationMultimap is an internal container.
namespace cppcrystal::data {

// rotation -> indices into operations_from_database(hall_number), built once
// per setting. Empty for an out-of-range Hall number.
[[nodiscard]] RotationMultimap<int> const &
operations_by_rotation(int hall_number);

} // namespace cppcrystal::data
