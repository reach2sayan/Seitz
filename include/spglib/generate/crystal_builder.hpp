#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/types.hpp>
#include <spglib/generate/distance_check.hpp>
#include <spglib/generate/wyckoff_combinations.hpp>
#include <spglib/group/point_group.hpp>
#include <spglib/group/space_group.hpp>

#include <cstdint>
#include <optional>

namespace spglib::generate {

// A generated crystal: the assembled conventional cell together with the Wyckoff
// assignment that produced it. The assignment's positions are non-owning
// pointers into `sg`, so the SpaceGroup passed to random_crystal must outlive
// the result (same contract as WyckoffCombination).
struct GeneratedCrystal {
  Cell cell;
  WyckoffCombination assignment;
};

// Controls for random_crystal. A struct rather than positional flags so callers
// name what they set (no magic argument order).
struct RandomCrystalOptions {
  // Multiplies the element-aware volume estimate (estimated_cell_volume).
  double volume_factor = 1.0;
  // Deterministic seed; std::nullopt selects a fixed default seed.
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / lattice resampling attempts per Wyckoff combination before
  // falling back to the next combination.
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion for a generated structure.
  DistanceTolerance distance = {};
  // Restrict placement to the general position (generic coordinates only). Useful
  // when the exact target group must be realized without the accidental extra
  // symmetry that special (fixed/high-symmetry) sites can introduce — notably for
  // layer groups, where atoms confined to a special site in one plane gain a
  // horizontal mirror. Requires the total atom count to be a multiple of the
  // general-position multiplicity.
  bool general_position_only = false;
};

// Generate a random crystal with the symmetry of `sg` and the given composition
// (PyXtal pyxtal.from_random / build). It enumerates the compatible Wyckoff
// assignments, samples them in a random order, and for each one repeatedly
// (`attempts_per_combination`) draws an element-aware random lattice and random
// free coordinates, accepting the first structure whose interatomic distances
// pass `options.distance`. Falls back to the next combination on exhaustion.
//
// Deterministic in `options.seed`. Errors with e_message when the composition is
// incompatible with the space group, or when no combination yields a
// distance-valid structure within the attempt budget.
[[nodiscard]] Result<GeneratedCrystal>
random_crystal(group::SpaceGroup const &sg, Composition const &comp,
               RandomCrystalOptions const &options = {});

// Generate a random 2D-periodic (layer) crystal with the symmetry of the layer
// group `lg` (a SpaceGroup built via SpaceGroup::from_layer_number / from_layer_
// hall) and the given composition. The structure is periodic in the a-b plane;
// the c axis is the non-periodic stacking direction, so the returned Cell has
// aperiodic_axis() == 2 (validate with spglib::get_layer_dataset). Atoms are
// placed on the layer Wyckoff positions, the in-plane cell is sized element-
// aware, and clashes are rejected with the aperiodic-aware distance check.
//
// Deterministic in options.seed. Errors when the composition is incompatible
// with the layer group, or when no distance-valid structure is found within the
// attempt budget.
[[nodiscard]] Result<GeneratedCrystal>
random_layer_crystal(group::SpaceGroup const &lg, Composition const &comp,
                     RandomCrystalOptions const &options = {});

// A generated 0D cluster: Cartesian atomic coordinates (row i = atom i) with
// their types, plus the metric `basis` that realises the point group's geometry.
// A point-group operation acts on the cluster in Cartesian space as
// `basis * rotation * basis.inverse()`; the cluster (as a point set) is invariant
// under every such operation. There is no lattice — the cluster is non-periodic.
struct GeneratedCluster {
  Positions coordinates;
  Types types;
  Matrix3d basis;
};

// Controls for random_cluster. A struct rather than positional flags so callers
// name what they set.
struct RandomClusterOptions {
  // Multiplies the element-aware size estimate (the cluster's characteristic
  // extent), so a larger factor spreads atoms further apart.
  double size_factor = 1.0;
  // Deterministic seed; std::nullopt selects a fixed default seed.
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / metric resampling attempts per Wyckoff combination before
  // falling back to the next combination.
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion (non-periodic, all-pairs).
  DistanceTolerance distance = {};
  // Restrict placement to the general position (generic points only), so the
  // cluster realises exactly the target point group with no accidental extra
  // symmetry from atoms confined to an axis or mirror plane. Requires the total
  // atom count to be a multiple of the general-position multiplicity.
  bool general_position_only = false;
};

// Generate a random 0D cluster with the symmetry of the point group `pg` and the
// given composition (PyXtal's molecular/cluster generation, dim = 0). It
// enumerates the compatible Wyckoff assignments, samples them in a random order,
// and for each one repeatedly draws a metric-compatible basis (so non-cubic
// point groups get the right geometry) and random free coordinates, accepting
// the first cluster whose interatomic distances pass `options.distance`.
//
// Deterministic in `options.seed`. Errors when the composition is incompatible
// with the point group's Wyckoff positions, or when no combination yields a
// distance-valid cluster within the attempt budget.
[[nodiscard]] Result<GeneratedCluster>
random_cluster(group::PointGroup const &pg, Composition const &comp,
               RandomClusterOptions const &options = {});

} // namespace spglib::generate
