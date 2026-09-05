#pragma once

#include <seitz/data/spg_database.hpp>

#include "core/matrix_order.hpp"

// The rotation index over a database setting's operations. Private: only the
// Hall-symbol matcher needs it, and RotationMultimap is an internal container.
namespace seitz::data {

// rotation -> indices into operations_from_database(hall), built once per
// setting.
[[nodiscard]] RotationMultimap<int> const &
operations_by_rotation(HallNumber hall);

} // namespace seitz::data
