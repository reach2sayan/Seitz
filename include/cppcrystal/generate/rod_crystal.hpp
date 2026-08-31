#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/group/rod_group.hpp>

namespace cppcrystal::generate {

// A generated 1D-periodic crystal; the cell is self-describing (its
// periodicity is {aperiodic, aperiodic, periodic} along the rod axis).
struct GeneratedRodCrystal {
  Cell cell;
};

// Generate a random 1D-periodic crystal with the symmetry of the rod group `rg`
// and the given composition. The structure is periodic along rg.periodic_axis()
// (c); the other two axes are aperiodic (a vacuum-padded rod cross-section).
// Atoms are placed on the rod Wyckoff positions and the orbit is expanded
// folding ONLY the periodic axis — a flipping operation along an aperiodic axis
// must land at the Cartesian image (e.g. -y), never the wrapped 1 - y (the rod
// analogue of the layer-group aperiodic-axis fix).
//
// Deterministic in options.seed. Errors when the composition is incompatible
// with the rod group's (currently general-position-only) Wyckoff positions, or
// when no distance-valid structure is found within the attempt budget.
[[nodiscard]] Result<GeneratedRodCrystal>
random_rod_crystal(group::RodGroup const &rg, Composition const &comp,
                   GenerateOptions const &options = {});

} // namespace cppcrystal::generate
