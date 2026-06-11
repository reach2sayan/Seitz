#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/generate/wyckoff_combinations.hpp>
#include <spglib/group/space_group.hpp>

#include <cstdint>
#include <optional>

namespace spglib::generate {

// Generate a random crystal structure with the symmetry of `sg` and the given
// composition (PyXtal pyxtal.from_random / build). Picks a compatible Wyckoff
// assignment, generates a crystal-system-consistent random lattice sized for the
// atom count (scaled by `volume_factor`), places each chosen Wyckoff position at
// a random free coordinate, and returns the assembled conventional Cell.
//
// Deterministic in `seed` (std::nullopt selects a fixed default seed). Errors
// with e_message when no Wyckoff assignment realizes the composition.
[[nodiscard]] Result<Cell>
random_crystal(group::SpaceGroup const &sg, Composition const &comp,
               double volume_factor = 1.0,
               std::optional<std::uint64_t> seed = std::nullopt);

} // namespace spglib::generate
