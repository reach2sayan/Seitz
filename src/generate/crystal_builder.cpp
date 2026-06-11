#include <spglib/generate/crystal_builder.hpp>

#include <spglib/generate/random_lattice.hpp>

#include <spglib/math/fractional.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace spglib::generate {

namespace {

// Expand one placement's orbit for a layer cell: the coset operations are
// applied without folding the non-periodic axis (folding it would send a
// c-flipping image to 1-z instead of -z and break that operation), and only the
// periodic axes are wrapped to [0, 1).
void expand_layer_orbit(WyckoffCombination::Placement const &placement,
                        Vector3d const &xyz, int aperiodic_axis,
                        std::vector<Vector3d> &rows, Types &types) {
  Vector3d const canonical = placement.position->canonical_coordinate(xyz);
  for (auto const &op : placement.position->orbit_operations()) {
    Vector3d image = op.apply(canonical);
    for (int axis = 0; axis < 3; ++axis) {
      if (axis != aperiodic_axis) {
        image[axis] = math::wrap_to_unit_cell(image[axis]);
      }
    }
    rows.push_back(image);
    types.push_back(placement.type);
  }
}

// Assemble the conventional cell for one Wyckoff assignment: each placement's
// free coordinate is drawn uniformly and expanded into its full orbit.
// `aperiodic_axis` is carried into the Cell (set for layer-group generation),
// where the orbit is built unwrapped along that axis and the layer is then
// centred in the cell.
[[nodiscard]] Cell assemble(WyckoffCombination const &combo,
                            Matrix3d const &lattice, std::mt19937_64 &rng,
                            std::optional<int> aperiodic_axis = std::nullopt) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placement : combo.placements) {
    if (aperiodic_axis) {
      // The free aperiodic coordinate is drawn near 0 (a thin band) so a
      // c-flipping image lands at -z, also near 0, keeping the layer thin; the
      // layer is recentred below.
      constexpr double kBand = 0.05;
      Vector3d xyz{unit(rng), unit(rng), unit(rng)};
      xyz[*aperiodic_axis] = kBand * (2.0 * unit(rng) - 1.0);
      expand_layer_orbit(placement, xyz, *aperiodic_axis, rows, types);
    } else {
      Vector3d const xyz{unit(rng), unit(rng), unit(rng)};
      Positions const orbit = placement.position->get_all_positions(xyz);
      for (Index i = 0; i < orbit.rows(); ++i) {
        rows.push_back(orbit.row(i).transpose());
        types.push_back(placement.type);
      }
    }
  }

  // The orbit is built symmetric about the operations' fixed planes (near c = 0,
  // from the small seed band); it is left there rather than recentred, which
  // would move the layer's symmetry planes away from the database operations.
  Positions positions(static_cast<Index>(rows.size()), 3);
  for (Index i = 0; i < static_cast<Index>(rows.size()); ++i) {
    positions.row(i) = rows[static_cast<std::size_t>(i)].transpose();
  }
  return Cell{lattice, std::move(positions), std::move(types), aperiodic_axis};
}

// Drop every combination that uses anything but the general position, when the
// caller asked for general-position-only generation.
void restrict_to_general_position(std::vector<WyckoffCombination> &combos,
                                  group::SpaceGroup const &group) {
  group::WyckoffPosition const *general = &group.wyckoffs().back();
  std::erase_if(combos, [&](WyckoffCombination const &combo) {
    return std::ranges::any_of(combo.placements, [&](auto const &placement) {
      return placement.position != general;
    });
  });
}

} // namespace

Result<GeneratedCrystal> random_crystal(group::SpaceGroup const &sg,
                                        Composition const &comp,
                                        RandomCrystalOptions const &options) {
  std::vector<WyckoffCombination> combos = list_wyckoff_combinations(sg, comp);
  if (combos.empty()) {
    return leaf::new_error(e_message{
        "generate::random_crystal: composition is not compatible with the "
        "Wyckoff positions of the requested space group"});
  }

  if (options.general_position_only) {
    restrict_to_general_position(combos, sg);
    if (combos.empty()) {
      return leaf::new_error(e_message{
          "generate::random_crystal: general_position_only requires the atom "
          "count to be a multiple of the general-position multiplicity"});
    }
  }

  std::mt19937_64 rng(options.seed.value_or(0));

  // Strategic sampling: try the compatible assignments in a random order rather
  // than always the first/one blindly-indexed combination.
  std::ranges::shuffle(combos, rng);

  double const target_volume =
      options.volume_factor * estimated_cell_volume(comp);
  CrystalSystem const system = crystal_system(sg.number());

  for (auto const &combo : combos) {
    for (int attempt = 0; attempt < options.attempts_per_combination;
         ++attempt) {
      Matrix3d const lattice = random_lattice(system, target_volume, rng());
      Cell cell = assemble(combo, lattice, rng);
      if (distances_valid(cell, options.distance)) {
        return GeneratedCrystal{std::move(cell), combo};
      }
    }
  }

  return leaf::new_error(e_message{
      "generate::random_crystal: no distance-valid structure found within the "
      "attempt budget for any compatible Wyckoff assignment"});
}

Result<GeneratedCrystal>
random_layer_crystal(group::SpaceGroup const &lg, Composition const &comp,
                     RandomCrystalOptions const &options) {
  std::vector<WyckoffCombination> combos = list_wyckoff_combinations(lg, comp);
  if (combos.empty()) {
    return leaf::new_error(e_message{
        "generate::random_layer_crystal: composition is not compatible with the "
        "Wyckoff positions of the requested layer group"});
  }

  if (options.general_position_only) {
    restrict_to_general_position(combos, lg);
    if (combos.empty()) {
      return leaf::new_error(e_message{
          "generate::random_layer_crystal: general_position_only requires the "
          "atom count to be a multiple of the general-position multiplicity"});
    }
  }

  std::mt19937_64 rng(options.seed.value_or(0));
  std::ranges::shuffle(combos, rng);

  // The layer is periodic in the a-b plane; c (aperiodic axis 2) is the stacking
  // direction with a large vacuum gap so the thin layer is isolated. The in-plane
  // area is sized element-aware (the atoms occupy a slab a few angstrom thick).
  constexpr double kLayerThickness = 3.0; // angstrom, for the in-plane area
  constexpr double kVacuum = 18.0;        // angstrom, the c (aperiodic) length
  double const area =
      options.volume_factor * estimated_cell_volume(comp) / kLayerThickness;
  double const c_length = kVacuum;

  for (auto const &combo : combos) {
    for (int attempt = 0; attempt < options.attempts_per_combination;
         ++attempt) {
      Matrix3d const lattice =
          random_layer_lattice(lg.operations(), area, c_length, rng());
      Cell cell = assemble(combo, lattice, rng, /*aperiodic_axis=*/2);
      if (distances_valid(cell, options.distance)) {
        return GeneratedCrystal{std::move(cell), combo};
      }
    }
  }

  return leaf::new_error(e_message{
      "generate::random_layer_crystal: no distance-valid structure found within "
      "the attempt budget for any compatible Wyckoff assignment"});
}

} // namespace spglib::generate
