#include <cppcrystal/generate/generator.hpp>

#include "generate/random_lattice.hpp"
#include <cppcrystal/data/element_data.hpp>
#include <cppcrystal/generate/distance_check.hpp>

#include <algorithm>
#include <cmath>
#include <format>
#include <ranges>
#include <utility>
#include <vector>

namespace cppcrystal::generate {

namespace {

// A layer sits in a slab a few angstrom thick, padded along c by a vacuum gap
// wide enough to isolate it from its own image; a rod is padded the same way
// across its two aperiodic axes.
constexpr double kLayerThickness = 3.0; // angstrom, for the in-plane area
constexpr double kVacuum = 18.0;        // angstrom, the aperiodic padding

// The 1D counterpart of estimated_cell_volume: a rod's atoms stack along its
// periodic axis, so the repeat must fit their summed diameters — a volume
// estimate would be nonsense here.
[[nodiscard]] double estimated_repeat(Composition const &comp,
                                      GenerateOptions const &options) {
  return std::ranges::fold_left(comp, 0.0, [&](double sum, auto const &entry) {
    auto const &[type, count] = entry;
    double const r =
        data::covalent_radius(type).value_or(options.distance.fallback_radius);
    return sum + static_cast<double>(count) * 2.0 * r;
  });
}

// One structure for one assignment: each placement's seed is drawn from the
// family's box, projected onto its Wyckoff locus and expanded into the full
// orbit, folding only the periodic axes.
[[nodiscard]] Cell assemble(Assignment<group::Wyckoff> const &assignment,
                            Matrix3d const &lattice,
                            CellPeriodicity const &periodicity,
                            SeedBox const &box, std::mt19937_64 &rng) {
  std::vector<Vector3d> rows;
  Types types;
  for (auto const &placed : assignment) {
    Positions const orbit =
        placed.position->orbit(box.sample(rng), periodicity);
    for (Index i = 0; i < orbit.rows(); ++i) {
      rows.push_back(orbit.row(i).transpose());
      types.push_back(placed.type);
    }
  }
  return Cell{Lattice{lattice}, to_positions(rows), std::move(types),
              periodicity};
}

// A seed box spanning the full repeat along every periodic axis and a band
// centred on 0, of half-width `band`, along every aperiodic one.
[[nodiscard]] SeedBox box_for(CellPeriodicity const &periodicity, double band) {
  SeedBox box;
  for (auto const [axis, kind] : periodicity | std::views::enumerate) {
    if (kind == AxisKind::aperiodic) {
      box.low[axis] = -band;
      box.high[axis] = band;
    }
  }
  return box;
}

} // namespace

Vector3d SeedBox::sample(std::mt19937_64 &rng) const {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  Vector3d const u{unit(rng), unit(rng), unit(rng)};
  return low + (high - low).cwiseProduct(u);
}

std::string_view
GroupTraits<group::SpaceGroup>::kind(group::SpaceGroup const &g) noexcept {
  return g.hall().family() == GroupFamily::layer ? "layer group"
                                                 : "space group";
}

CellPeriodicity GroupTraits<group::SpaceGroup>::periodicity(
    group::SpaceGroup const &g) noexcept {
  // A layer group is periodic in the a-b plane; c is the stacking direction.
  return g.hall().family() == GroupFamily::layer ? aperiodic_along(2)
                                                 : all_periodic();
}

SeedBox
GroupTraits<group::SpaceGroup>::seed_box(group::SpaceGroup const &g) noexcept {
  // The free aperiodic coordinate is drawn from a thin band so a c-flipping
  // image lands near 0 too and the layer stays thin, symmetric about the
  // operations' fixed planes. It is left there rather than recentred, which
  // would move those planes away from the database operations.
  return box_for(periodicity(g), 0.05);
}

Matrix3d GroupTraits<group::SpaceGroup>::lattice(group::SpaceGroup const &g,
                                                 Composition const &comp,
                                                 GenerateOptions const &options,
                                                 int /*attempt*/,
                                                 std::mt19937_64 &rng) {
  double const volume = options.scale * estimated_cell_volume(comp);
  if (g.hall().family() == GroupFamily::layer) {
    // random_layer_lattice symmetrizes the in-plane metric over the operations'
    // 2x2 blocks, so it is exactly invariant under the in-plane point group
    // whatever the crystal system — no number-range table, which would
    // mis-handle mixed cases (p112 wants an oblique gamma where pm11 forces
    // gamma = 90 in the same range).
    return random_layer_lattice(g.operations(), volume / kLayerThickness,
                                kVacuum, rng());
  }
  return random_lattice(crystal_system(g.number()), volume, rng());
}

std::string_view
GroupTraits<group::PointGroup>::kind(group::PointGroup const &) noexcept {
  return "point group";
}

CellPeriodicity GroupTraits<group::PointGroup>::periodicity(
    group::PointGroup const &) noexcept {
  return none_periodic();
}

SeedBox
GroupTraits<group::PointGroup>::seed_box(group::PointGroup const &) noexcept {
  return box_for(none_periodic(),
                 0.5); // the whole metric, centred on the origin
}

Matrix3d GroupTraits<group::PointGroup>::lattice(group::PointGroup const &g,
                                                 Composition const &comp,
                                                 GenerateOptions const &options,
                                                 int attempt,
                                                 std::mt19937_64 &rng) {
  // The point group's crystal system fixes the metric the cluster is realized
  // in, so trigonal/hexagonal clusters get the correct non-orthogonal geometry
  // and the integer operations become true Cartesian isometries. The metric
  // grows on repeated failure, to spread the atoms apart.
  double const volume = options.scale * estimated_cell_volume(comp) *
                        (1.0 + 0.2 * static_cast<double>(attempt));
  return random_lattice(crystal_system(g.representative_spacegroup()), volume,
                        rng());
}

std::string_view
GroupTraits<group::RodGroup>::kind(group::RodGroup const &) noexcept {
  return "rod group";
}

CellPeriodicity
GroupTraits<group::RodGroup>::periodicity(group::RodGroup const &g) noexcept {
  return periodic_along(g.periodic_axis());
}

SeedBox
GroupTraits<group::RodGroup>::seed_box(group::RodGroup const &g) noexcept {
  // The cross-section band is centred on 0 (so an axis-flipping operation maps
  // onto a real atom and the rod stays centred) but NOT thin: a too-narrow band
  // collapses the images of a rotation about the periodic axis onto the axis
  // itself, where they clash.
  return box_for(periodicity(g), 0.25);
}

Matrix3d GroupTraits<group::RodGroup>::lattice(group::RodGroup const &g,
                                               Composition const &comp,
                                               GenerateOptions const &options,
                                               int attempt,
                                               std::mt19937_64 &rng) {
  // random_layer_lattice symmetrizes the (a, b) metric and sets the c length:
  // exactly a rod's geometry, only the periodic/aperiodic reading differs.
  double const repeat = options.scale * estimated_repeat(comp, options) *
                        (1.0 + 0.2 * static_cast<double>(attempt));
  return random_layer_lattice(g.operations(), kVacuum * kVacuum, repeat, rng());
}

template <class G>
Result<Generated> Generator<G>::operator()(Composition const &comp) const {
  using Traits = GroupTraits<G>;
  auto const positions = group_->wyckoffs();
  auto const kind = Traits::kind(*group_);

  // Not const: shuffled below.
  std::vector<Assignment<group::Wyckoff>> assignments(
      std::from_range, enumerate_assignments(positions, comp) |
                           std::views::take(kMaxAssignments));
  if (assignments.empty()) {
    return leaf::new_error(e_message{
        std::format("generate: the composition is not compatible with the "
                    "Wyckoff positions of the requested {}",
                    kind)});
  }

  if (options_.placement == Placement::general_only) {
    group::Wyckoff const *general = &positions.back();
    std::erase_if(assignments, [&](Assignment<group::Wyckoff> const &a) {
      return std::ranges::any_of(a, [&](Placed<group::Wyckoff> const &p) {
        return p.position != general;
      });
    });
    if (assignments.empty()) {
      return leaf::new_error(e_message{std::format(
          "generate: Placement::general_only requires the atom count to be a "
          "multiple of the general-position multiplicity of the requested {}",
          kind)});
    }
  }

  std::mt19937_64 rng(options_.seed.value_or(0));
  // Strategic sampling: try the compatible assignments in a random order rather
  // than always the first/one blindly-indexed combination.
  std::ranges::shuffle(assignments, rng);

  CellPeriodicity const periodicity = Traits::periodicity(*group_);
  SeedBox const box = Traits::seed_box(*group_);
  for (auto const &assignment : assignments) {
    for (int attempt = 0; attempt < options_.attempts_per_combination;
         ++attempt) {
      Matrix3d const lattice =
          Traits::lattice(*group_, comp, options_, attempt, rng);
      Cell cell = assemble(assignment, lattice, periodicity, box, rng);
      if (distances_valid(cell, options_.distance)) {
        return Generated{std::move(cell), assignment};
      }
    }
  }

  return leaf::new_error(e_message{std::format(
      "generate: no distance-valid structure found within the attempt budget "
      "for any compatible Wyckoff assignment on the requested {}",
      kind)});
}

template class Generator<group::SpaceGroup>;
template class Generator<group::PointGroup>;
template class Generator<group::RodGroup>;

} // namespace cppcrystal::generate
