#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/periodicity.hpp>
#include <spglib/generate/distance_check.hpp>
#include <spglib/generate/wyckoff_combinations.hpp>
#include <spglib/group/rod_group.hpp>

#include <cstdint>
#include <optional>

namespace spglib::generate {

struct GeneratedRodCrystal {
  Cell cell;
  CellPeriodicity periodicity;
};

// Controls for random_rod_crystal (named fields, no positional flags).
struct RandomRodOptions {
  // Multiplies the element-aware size estimate (sets the periodic repeat length
  // and the vacuum cross-section).
  double size_factor = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  int attempts_per_combination = 50;
  DistanceTolerance distance = {};
  bool general_position_only = false;
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
                   RandomRodOptions const &options = {});

} // namespace spglib::generate
