#include <cppcrystal/generate/crystal_builder.hpp>

#include <cppcrystal/generate/random_lattice.hpp>

#include "math/fractional.hpp"

#include <array>
#include <optional>
#include <random>
#include <span>
#include <vector>

namespace cppcrystal::generate {

namespace {

// Assemble the conventional cell for one Wyckoff assignment: each placement's
// free coordinate is drawn uniformly and expanded into its full orbit.
// `aperiodic_axis` is carried into the Cell (set for layer-group generation),
// where the orbit is built unwrapped along that axis and the layer is then
// centred in the cell.
[[nodiscard]] Cell assemble(WyckoffCombination const &combo,
                            Matrix3d const &lattice, std::mt19937_64 &rng,
                            std::optional<int> aperiodic_axis = std::nullopt) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  CellPeriodicity const periodicity = aperiodic_axis
                                          ? aperiodic_along(*aperiodic_axis)
                                          : all_periodic();
  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placement : combo) {
    Vector3d xyz{unit(rng), unit(rng), unit(rng)};
    if (aperiodic_axis) {
      // The free aperiodic coordinate is drawn near 0 (a thin band) so a
      // c-flipping image lands at -z, also near 0, keeping the layer thin; the
      // layer is recentred below. Wyckoff::orbit leaves that axis unfolded.
      constexpr double kBand = 0.05;
      xyz[*aperiodic_axis] = kBand * (2.0 * unit(rng) - 1.0);
    }
    Positions const orbit = placement.position->orbit(xyz, periodicity);
    for (Index i = 0; i < orbit.rows(); ++i) {
      rows.push_back(orbit.row(i).transpose());
      types.push_back(placement.type);
    }
  }

  // The orbit is built symmetric about the operations' fixed planes (near c =
  // 0, from the small seed band); it is left there rather than recentred, which
  // would move the layer's symmetry planes away from the database operations.
  return Cell{Lattice{lattice}, to_positions(rows), std::move(types),
              aperiodic_axis ? aperiodic_along(*aperiodic_axis)
                             : all_periodic()};
}

} // namespace

Result<GeneratedCrystal> random_crystal(group::SpaceGroup const &sg,
                                        Composition const &comp,
                                        GenerateOptions const &options) {
  double const target_volume = options.scale * estimated_cell_volume(comp);
  CrystalSystem const system = crystal_system(sg.number());

  return detail::search_assignments(
      sg.wyckoffs(), comp, options, "generate::random_crystal", "space group",
      [&](WyckoffCombination const &combo, int /*attempt*/,
          std::mt19937_64 &rng) -> std::optional<GeneratedCrystal> {
        Matrix3d const lattice = random_lattice(system, target_volume, rng());
        Cell cell = assemble(combo, lattice, rng);
        if (!distances_valid(cell, options.distance)) {
          return std::nullopt;
        }
        return GeneratedCrystal{std::move(cell), combo};
      });
}

Result<GeneratedCrystal>
random_layer_crystal(group::SpaceGroup const &lg, Composition const &comp,
                     GenerateOptions const &options) {
  // The layer is periodic in the a-b plane; c (aperiodic axis 2) is the
  // stacking direction with a large vacuum gap so the thin layer is isolated.
  // The in-plane area is sized element-aware (the atoms occupy a slab a few
  // angstrom thick).
  constexpr double kLayerThickness = 3.0; // angstrom, for the in-plane area
  constexpr double kVacuum = 18.0;        // angstrom, the c (aperiodic) length
  double const area =
      options.scale * estimated_cell_volume(comp) / kLayerThickness;
  double const c_length = kVacuum;

  return detail::search_assignments(
      lg.wyckoffs(), comp, options, "generate::random_layer_crystal",
      "layer group",
      [&](WyckoffCombination const &combo, int /*attempt*/,
          std::mt19937_64 &rng) -> std::optional<GeneratedCrystal> {
        Matrix3d const lattice =
            random_layer_lattice(lg.operations(), area, c_length, rng());
        Cell cell = assemble(combo, lattice, rng, /*aperiodic_axis=*/2);
        if (!distances_valid(cell, options.distance)) {
          return std::nullopt;
        }
        return GeneratedCrystal{std::move(cell), combo};
      });
}

namespace {

// Assemble one cluster: each placement's free coordinate is drawn uniformly in
// fractional space and expanded over its orbit, then mapped to Cartesian by the
// metric basis. The result is non-periodic.
[[nodiscard]] std::pair<Positions, Types>
assemble_cluster(Assignment<group::Wyckoff> const &combo,
                 Matrix3d const &basis, std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> coord(-0.5, 0.5);

  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placement : combo) {
    int const dof = placement.position->degrees_of_freedom();
    std::array<double, 3> params{};
    for (int i = 0; i < dof; ++i) {
      params[static_cast<std::size_t>(i)] = coord(rng);
    }
    Vector3d const q = placement.position->sample(
        std::span<double const>(params.data(), static_cast<std::size_t>(dof)));
    for (auto const &op : placement.position->orbit_operations()) {
      Vector3d const frac = op.rotation.cast<double>() * q;
      rows.push_back(basis * frac);
      types.push_back(placement.type);
    }
  }

  return {to_positions(rows), std::move(types)};
}

} // namespace

Result<GeneratedCluster> random_cluster(group::PointGroup const &pg,
                                        Composition const &comp,
                                        GenerateOptions const &options) {
  // The point group's crystal system fixes which metric the cluster is
  // realised in (so trigonal/hexagonal clusters get the correct,
  // non-orthogonal geometry and the integer operations become true Cartesian
  // isometries).
  CrystalSystem const system = crystal_system(pg.representative_spacegroup());
  double const base_volume = options.scale * estimated_cell_volume(comp);

  return detail::search_assignments(
      pg.wyckoffs(), comp, options, "generate::random_cluster", "point group",
      [&](Assignment<group::Wyckoff> const &combo, int attempt,
          std::mt19937_64 &rng) -> std::optional<GeneratedCluster> {
        // Grow the metric on repeated failure to spread atoms apart.
        double const volume = base_volume * (1.0 + 0.2 * attempt);
        Matrix3d const basis = random_lattice(system, volume, rng());
        auto [coordinates, types] = assemble_cluster(combo, basis, rng);
        if (!cluster_distances_valid(coordinates, types, options.distance)) {
          return std::nullopt;
        }
        return GeneratedCluster{std::move(coordinates), std::move(types),
                                basis};
      });
}

} // namespace cppcrystal::generate
