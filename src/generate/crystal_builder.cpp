#include <spglib/generate/crystal_builder.hpp>

#include <spglib/generate/random_lattice.hpp>

#include <random>
#include <vector>

namespace spglib::generate {

namespace {
// Typical atomic volume (cubic angstrom per atom) used to size the cell before
// the caller's volume_factor is applied.
constexpr double kVolumePerAtom = 20.0;
} // namespace

Result<Cell> random_crystal(group::SpaceGroup const &sg,
                            Composition const &comp, double volume_factor,
                            std::optional<std::uint64_t> seed) {
  std::vector<WyckoffCombination> const combos =
      list_wyckoff_combinations(sg, comp);
  if (combos.empty()) {
    return leaf::new_error(e_message{
        "generate::random_crystal: composition is not compatible with the "
        "Wyckoff positions of the requested space group"});
  }

  std::mt19937_64 rng(seed.value_or(0));
  std::uniform_real_distribution<double> unit(0.0, 1.0);

  WyckoffCombination const &combo = combos[rng() % combos.size()];

  int total_atoms = 0;
  for (auto const &[type, count] : comp) {
    total_atoms += count;
  }
  double const target_volume =
      volume_factor * static_cast<double>(total_atoms) * kVolumePerAtom;
  Matrix3d const lattice =
      random_lattice(crystal_system(sg.number()), target_volume, rng());

  std::vector<Vector3d> rows;
  Types types;
  rows.reserve(static_cast<std::size_t>(total_atoms));
  types.reserve(static_cast<std::size_t>(total_atoms));
  for (auto const &placement : combo.placements) {
    Vector3d const xyz{unit(rng), unit(rng), unit(rng)};
    Positions const orbit = placement.position->get_all_positions(xyz);
    for (Index i = 0; i < orbit.rows(); ++i) {
      rows.push_back(orbit.row(i).transpose());
      types.push_back(placement.type);
    }
  }

  Positions positions(static_cast<Index>(rows.size()), 3);
  for (Index i = 0; i < static_cast<Index>(rows.size()); ++i) {
    positions.row(i) = rows[static_cast<std::size_t>(i)].transpose();
  }

  return Cell{lattice, std::move(positions), std::move(types)};
}

} // namespace spglib::generate
