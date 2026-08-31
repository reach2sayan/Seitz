#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/generate/assignments.hpp>
#include <cppcrystal/generate/wyckoff_combinations.hpp>
#include <cppcrystal/group/point_group.hpp>
#include <cppcrystal/group/space_group.hpp>

namespace cppcrystal::generate {

struct GeneratedCrystal {
  Cell cell;
  WyckoffCombination assignment;
};

[[nodiscard]] Result<GeneratedCrystal>
random_crystal(group::SpaceGroup const &sg, Composition const &comp,
               GenerateOptions const &options = {});

[[nodiscard]] Result<GeneratedCrystal>
random_layer_crystal(group::SpaceGroup const &lg, Composition const &comp,
                     GenerateOptions const &options = {});

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

[[nodiscard]] Result<GeneratedCluster>
random_cluster(group::PointGroup const &pg, Composition const &comp,
               GenerateOptions const &options = {});

} // namespace cppcrystal::generate
