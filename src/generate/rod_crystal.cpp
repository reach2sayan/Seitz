#include <cppcrystal/generate/rod_crystal.hpp>

#include <cppcrystal/data/element_data.hpp>
#include <cppcrystal/generate/random_lattice.hpp>
#include <cppcrystal/math/fractional.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <random>
#include <span>
#include <vector>

namespace cppcrystal::generate {

namespace {

// Assemble one rod cell: each placement's free coordinate is drawn with the
// periodic axis spanning [0, 1) and the aperiodic axes confined to a thin band
// near 0 (so a flipping operation's image stays near the rod core), then the
// orbit is expanded folding ONLY the periodic axis.
[[nodiscard]] Cell assemble_rod(Assignment<group::LocusWyckoff> const &combo,
                                Matrix3d const &lattice, int periodic_axis,
                                std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  // Aperiodic axes are sampled in a band CENTRED on 0 (so an axis-flipping op,
  // a -> -a, maps onto a real atom and the rod stays centred on the periodic
  // axis) but NOT thin: a too-narrow band collapses the images of a rotation
  // about the periodic axis onto the axis itself, where they clash. The
  // periodic axis spans the full [0, 1) repeat. (This is the rod counterpart
  // of the layer generator's thin-c band, with the roles of periodic/aperiodic
  // swapped.)
  constexpr double kCrossSection = 0.25;

  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placement : combo) {
    int const dof = placement.position->degrees_of_freedom();
    // Each free coordinate is sampled by the kind of its locus direction (the
    // basis is axis-separated): a direction along the periodic axis spans the
    // full [0, 1) repeat; an aperiodic (cross-section) direction is confined
    // to the centred band. The locus origin supplies any fixed offset (special
    // positions), so this works for the general position and the specials
    // alike.
    std::array<double, 3> params{};
    for (int i = 0; i < dof; ++i) {
      Vector3d const dir = placement.position->free_direction(i);
      bool const along_periodic = std::abs(dir[periodic_axis]) > 0.5;
      params[static_cast<std::size_t>(i)] =
          along_periodic ? unit(rng) : kCrossSection * (2.0 * unit(rng) - 1.0);
    }
    Vector3d const q = placement.position->sample(
        std::span<double const>(params.data(), static_cast<std::size_t>(dof)));

    for (auto const &op : placement.position->orbit_operations()) {
      Vector3d image = op.apply(q);
      image[periodic_axis] = math::wrap_to_unit_cell(image[periodic_axis]);
      rows.push_back(image);
      types.push_back(placement.type);
    }
  }

  return Cell{lattice, to_positions(rows), std::move(types)};
}

} // namespace

Result<GeneratedRodCrystal>
random_rod_crystal(group::RodGroup const &rg, Composition const &comp,
                   GenerateOptions const &options) {
  int const periodic_axis = rg.periodic_axis();

  // The rod is periodic along c with a vacuum-padded cross-section. The
  // periodic repeat is sized by a 1D (linear) estimate — the atoms stack along
  // c, so the repeat must fit their summed diameters, NOT a 3D volume — while
  // a large cross-section area isolates neighbouring rods.
  // random_layer_lattice symmetrises the cross-section (a, b) metric over the
  // operations' 2x2 blocks and sets the c (here periodic) length: exactly the
  // geometry a rod needs (only the periodic/aperiodic interpretation differs).
  constexpr double kVacuum = 18.0; // angstrom, the aperiodic padding scale
  double const linear = std::ranges::fold_left(
      comp, 0.0, [&](double sum, auto const &entry) {
        auto const &[type, count] = entry;
        double const r = data::covalent_radius(type).value_or(
            options.distance.fallback_radius);
        return sum + static_cast<double>(count) * 2.0 * r;
      });
  double const repeat = options.scale * linear; // 1D packing along c
  double const area = kVacuum * kVacuum;        // vacuum cross-section

  CellPeriodicity periodicity = all_periodic();
  for (int axis = 0; axis < 3; ++axis) {
    if (axis != periodic_axis) {
      periodicity[static_cast<std::size_t>(axis)] = AxisKind::aperiodic;
    }
  }

  return detail::search_assignments(
      rg.wyckoffs(), comp, options, "generate::random_rod_crystal",
      "rod group",
      [&](Assignment<group::LocusWyckoff> const &combo, int attempt,
          std::mt19937_64 &rng) -> std::optional<GeneratedRodCrystal> {
        double const c_length = repeat * (1.0 + 0.2 * attempt);
        Matrix3d const lattice =
            random_layer_lattice(rg.operations(), area, c_length, rng());
        Cell cell = assemble_rod(combo, lattice, periodic_axis, rng);
        // Self-describing: {aperiodic, aperiodic, periodic}.
        cell.set_periodicity(periodicity);
        if (!distances_valid(cell, options.distance)) {
          return std::nullopt;
        }
        return GeneratedRodCrystal{std::move(cell)};
      });
}

} // namespace cppcrystal::generate
