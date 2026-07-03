#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/generate/distance_check.hpp>
#include <cppcrystal/generate/wyckoff_combinations.hpp>
#include <cppcrystal/group/point_group.hpp>
#include <cppcrystal/group/space_group.hpp>

#include <cstdint>
#include <optional>

namespace cppcrystal::generate {

struct GeneratedCrystal {
  Cell cell;
  WyckoffCombination assignment;
};

struct RandomCrystalOptions {
  double volume_factor = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  // Free-coordinate / lattice resampling attempts per Wyckoff combination
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion for a generated structure.
  DistanceTolerance distance = {};
  // Restrict placement to the general position (generic coordinates only).
  // Useful when the exact target group must be realized without the accidental
  // extra symmetry that special (fixed/high-symmetry) sites can introduce —
  // notably for layer groups, where atoms confined to a special site in one
  // plane gain a horizontal mirror. Requires the total atom count to be a
  // multiple of the general-position multiplicity.
  bool general_position_only = false;
};

[[nodiscard]] Result<GeneratedCrystal>
random_crystal(group::SpaceGroup const &sg, Composition const &comp,
               RandomCrystalOptions const &options = {});

[[nodiscard]] Result<GeneratedCrystal>
random_layer_crystal(group::SpaceGroup const &lg, Composition const &comp,
                     RandomCrystalOptions const &options = {});

// A generated 0D cluster: Cartesian atomic coordinates (row i = atom i) with
// their types, plus the metric `basis` that realizes the point group's
// geometry. A point-group operation acts on the cluster in Cartesian space as
// `basis * rotation * basis.inverse()`; the cluster (as a point set) is
// invariant under every such operation. There is no lattice — the cluster is
// non-periodic.
struct GeneratedCluster {
  Positions coordinates;
  Types types;
  Matrix3d basis;
};

struct RandomClusterOptions {
  double size_factor = 1.0;
  std::optional<std::uint64_t> seed = std::nullopt;
  int attempts_per_combination = 50;
  // Minimum-distance acceptance criterion (non-periodic, all-pairs).
  DistanceTolerance distance = {};
  bool general_position_only = false;
};

[[nodiscard]] Result<GeneratedCluster>
random_cluster(group::PointGroup const &pg, Composition const &comp,
               RandomClusterOptions const &options = {});

} // namespace cppcrystal::generate
