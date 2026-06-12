#pragma once

#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>
#include <spglib/data/spg_database.hpp>

#include <optional>

namespace spglib::spacegroup {

// Try to match a set of symmetry operations (given in the bravais /
// conventional setting) against a specific Hall setting. Returns the origin
// shift (in the bravais setting) if the operations reproduce that Hall
// setting's database operations, std::nullopt otherwise. 3D space groups only.

[[nodiscard]] std::optional<Vector3d>
match_hall_symbol(Matrix3d const &bravais_lattice, int hall_number,
                  data::Centering centering, SymmetryOperations const &symmetry,
                  double symprec);

} // namespace spglib::spacegroup
