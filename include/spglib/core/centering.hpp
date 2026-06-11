#pragma once

#include <spglib/core/types.hpp>
#include <spglib/data/spg_database.hpp> // data::Centering

#include <vector>

// Centering change-of-basis data, the single source of truth shared by the
// space-group search (spacegroup.cpp / hall_symbol.cpp) and cell
// standardization (standardize_cell). For a centering C:
//   - centering_matrix(C)      M     : primitive -> conventional basis
//   (integer)
//   - centering_matrix_inv(C)  M^-1  : conventional -> primitive basis
//   - centering_shifts(C)            : the non-trivial centering translations
//                                      (multiplicity - 1 of them; the zero
//                                      shift is implicit)
// PRIMITIVE (and any unused entry) maps to the identity / an empty shift list.
namespace spglib {

[[nodiscard]] Matrix3i const &centering_matrix(data::Centering c);
[[nodiscard]] Matrix3d const &centering_matrix_inv(data::Centering c);
[[nodiscard]] std::vector<Vector3d> centering_shifts(data::Centering c);

} // namespace spglib
