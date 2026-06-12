#include <spglib/generate/rod_crystal.hpp>

#include <spglib/data/element_data.hpp>
#include <spglib/generate/random_lattice.hpp>
#include <spglib/math/fractional.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace spglib::generate {

namespace {

// --- Distance check, rod flavour --------------------------------------------
// Only the periodic axis is periodic (the rod inverse of distance_check, where
// at most one axis is aperiodic). Folds just that component; the two aperiodic
// axes keep their raw difference and are not searched over neighbouring cells.

[[nodiscard]] double rod_min_distance(Vector3d const &a, Vector3d const &b,
                                      Matrix3d const &lattice, int periodic_axis,
                                      bool include_origin) {
  Vector3d const diff = a - b;
  Vector3d base = diff;
  double const p = diff[periodic_axis];
  base[periodic_axis] = p - std::round(p); // fold into [-0.5, 0.5]

  double best = std::numeric_limits<double>::infinity();
  for (int n = -1; n <= 1; ++n) {
    Vector3d off = base;
    off[periodic_axis] = base[periodic_axis] + static_cast<double>(n);
    if (!include_origin && off.squaredNorm() < 1e-20) {
      continue;
    }
    double const d2 = (lattice * off).squaredNorm();
    if (d2 < best) {
      best = d2;
    }
  }
  return std::sqrt(best);
}

[[nodiscard]] bool rod_distances_valid(Cell const &cell, int periodic_axis,
                                       DistanceTolerance tol) {
  Index const n = cell.size();
  std::vector<double> radius(static_cast<std::size_t>(n));
  std::ranges::transform(cell.types(), radius.begin(), [&](int type) {
    return data::covalent_radius(type).value_or(tol.fallback_radius);
  });

  for (Index i = 0; i < n; ++i) {
    Vector3d const pi = cell.position(i);
    double const ri = radius[static_cast<std::size_t>(i)];
    for (Index j = i; j < n; ++j) {
      double const min_dist =
          tol.scale * (ri + radius[static_cast<std::size_t>(j)]);
      double const dist = rod_min_distance(pi, cell.position(j), cell.lattice(),
                                           periodic_axis, /*include_origin=*/i != j);
      if (dist < min_dist) {
        return false;
      }
    }
  }
  return true;
}

// --- Wyckoff-assignment enumeration (the rod analogue of the 0D/2D ones) -----

struct RodPlacement {
  int type;
  group::RodWyckoff const *position{};
};
using RodCombination = std::vector<RodPlacement>;

struct RodEnumerator {
  std::span<group::RodWyckoff const> positions;
  std::vector<std::pair<int, int>> const &elements; // (type, count)
  std::size_t max_combinations;

  std::vector<bool> used_special;
  RodCombination placements;
  std::vector<RodCombination> out;

  [[nodiscard]] bool full() const { return out.size() >= max_combinations; }

  void recurse(std::size_t elem, std::size_t pos, int remaining) {
    if (full()) {
      return;
    }
    if (remaining == 0) {
      if (elem + 1 == elements.size()) {
        out.push_back(placements);
      } else {
        recurse(elem + 1, 0, elements[elem + 1].second);
      }
      return;
    }
    if (pos >= positions.size()) {
      return;
    }
    group::RodWyckoff const &wp = positions[pos];
    int const mult = wp.multiplicity();
    bool const fixed = wp.degrees_of_freedom() == 0;
    int const max_copies =
        fixed ? (used_special[pos] ? 0 : 1) : remaining / mult;
    for (int copies = 0; copies <= max_copies && copies * mult <= remaining;
         ++copies) {
      if (full()) {
        return;
      }
      for (int k = 0; k < copies; ++k) {
        placements.push_back({elements[elem].first, &wp});
      }
      if (fixed && copies == 1) {
        used_special[pos] = true;
      }
      recurse(elem, pos + 1, remaining - copies * mult);
      if (fixed && copies == 1) {
        used_special[pos] = false;
      }
      for (int k = 0; k < copies; ++k) {
        placements.pop_back();
      }
    }
  }
};

[[nodiscard]] std::vector<RodCombination>
list_rod_combinations(group::RodGroup const &rg, Composition const &comp,
                      std::size_t max_combinations = 1000) {
  std::vector<std::pair<int, int>> elements;
  for (auto const &[type, count] : comp) {
    if (count > 0) {
      elements.emplace_back(type, count);
    }
  }
  if (elements.empty() || max_combinations == 0) {
    return {};
  }
  RodEnumerator e{rg.wyckoffs(),
                  elements,
                  max_combinations,
                  std::vector<bool>(rg.wyckoffs().size(), false),
                  {},
                  {}};
  e.recurse(0, 0, elements.front().second);
  return std::move(e.out);
}

// Assemble one rod cell: each placement's free coordinate is drawn with the
// periodic axis spanning [0, 1) and the aperiodic axes confined to a thin band
// near 0 (so a flipping operation's image stays near the rod core), then the
// orbit is expanded folding ONLY the periodic axis.
[[nodiscard]] Cell assemble_rod(RodCombination const &combo,
                                Matrix3d const &lattice, int periodic_axis,
                                std::mt19937_64 &rng) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  // Aperiodic axes are sampled in a band CENTRED on 0 (so an axis-flipping op,
  // a -> -a, maps onto a real atom and the rod stays centred on the periodic
  // axis) but NOT thin: a too-narrow band collapses the images of a rotation
  // about the periodic axis onto the axis itself, where they clash. The periodic
  // axis spans the full [0, 1) repeat. (This is the rod counterpart of the layer
  // generator's thin-c band, with the roles of periodic/aperiodic swapped.)
  constexpr double kCrossSection = 0.25;

  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placement : combo) {
    int const dof = placement.position->degrees_of_freedom();
    // SKELETON: only the general position (dof 3, identity locus) exists, so the
    // sampled parameters map directly to (x, y, z). When special positions land,
    // build the seed on the locus and respect the periodic axis there too.
    std::array<double, 3> params{};
    for (int i = 0; i < dof; ++i) {
      params[static_cast<std::size_t>(i)] =
          (i == periodic_axis) ? unit(rng)
                               : kCrossSection * (2.0 * unit(rng) - 1.0);
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

  Positions positions(static_cast<Index>(rows.size()), 3);
  for (Index i = 0; i < static_cast<Index>(rows.size()); ++i) {
    positions.row(i) = rows[static_cast<std::size_t>(i)].transpose();
  }
  return Cell{lattice, std::move(positions), std::move(types)};
}

} // namespace

Result<GeneratedRodCrystal>
random_rod_crystal(group::RodGroup const &rg, Composition const &comp,
                   RandomRodOptions const &options) {
  std::vector<RodCombination> combos = list_rod_combinations(rg, comp);
  if (combos.empty()) {
    return leaf::new_error(e_message{
        "generate::random_rod_crystal: composition is not compatible with the "
        "Wyckoff positions of the requested rod group"});
  }

  // With only the general position available, every combination already uses it;
  // the flag is honoured trivially. (Kept for parity / future special positions.)
  if (options.general_position_only) {
    group::RodWyckoff const *general = &rg.wyckoffs().back();
    std::erase_if(combos, [&](RodCombination const &combo) {
      return std::ranges::any_of(combo, [&](RodPlacement const &p) {
        return p.position != general;
      });
    });
    if (combos.empty()) {
      return leaf::new_error(e_message{
          "generate::random_rod_crystal: general_position_only requires the atom "
          "count to be a multiple of the general-position multiplicity"});
    }
  }

  std::mt19937_64 rng(options.seed.value_or(0));
  std::ranges::shuffle(combos, rng);

  int const periodic_axis = rg.periodic_axis();

  // The rod is periodic along c with a vacuum-padded cross-section. The periodic
  // repeat is sized by a 1D (linear) estimate — the atoms stack along c, so the
  // repeat must fit their summed diameters, NOT a 3D volume — while a large
  // cross-section area isolates neighbouring rods. random_layer_lattice
  // symmetrises the cross-section (a, b) metric over the operations' 2x2 blocks
  // and sets the c (here periodic) length: exactly the geometry a rod needs
  // (only the periodic/aperiodic interpretation differs).
  constexpr double kVacuum = 18.0; // angstrom, the aperiodic padding scale
  double linear = 0.0;
  for (auto const &[type, count] : comp) {
    double const r =
        data::covalent_radius(type).value_or(options.distance.fallback_radius);
    linear += static_cast<double>(count) * 2.0 * r;
  }
  double const repeat = options.size_factor * linear; // 1D packing along c
  double const area = kVacuum * kVacuum;               // vacuum cross-section

  CellPeriodicity periodicity = all_periodic();
  for (int axis = 0; axis < 3; ++axis) {
    if (axis != periodic_axis) {
      periodicity[static_cast<std::size_t>(axis)] = AxisKind::aperiodic;
    }
  }

  for (auto const &combo : combos) {
    for (int attempt = 0; attempt < options.attempts_per_combination;
         ++attempt) {
      double const c_length = repeat * (1.0 + 0.2 * attempt);
      Matrix3d const lattice =
          random_layer_lattice(rg.operations(), area, c_length, rng());
      Cell cell = assemble_rod(combo, lattice, periodic_axis, rng);
      if (rod_distances_valid(cell, periodic_axis, options.distance)) {
        return GeneratedRodCrystal{std::move(cell), periodicity};
      }
    }
  }

  return leaf::new_error(e_message{
      "generate::random_rod_crystal: no distance-valid structure found within "
      "the attempt budget for any compatible Wyckoff assignment"});
}

} // namespace spglib::generate
