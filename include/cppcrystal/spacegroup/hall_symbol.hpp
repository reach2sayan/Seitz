#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/spg_database.hpp>

#include <optional>

namespace cppcrystal::spacegroup {

// Try to match a set of symmetry operations (given in the bravais /
// conventional setting) against a specific Hall setting. Returns the origin
// shift (in the bravais setting) if the operations reproduce that Hall
// setting's database operations, std::nullopt otherwise. 3D space groups only.

[[nodiscard]] std::optional<Vector3d>
match_hall_symbol(Matrix3d const &bravais_lattice, int hall_number,
                  data::Centering centering, SymmetryOperations const &symmetry,
                  double symprec);

} // namespace cppcrystal::spacegroup
