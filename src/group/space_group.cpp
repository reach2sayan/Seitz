#include <cppcrystal/group/space_group.hpp>

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/group/detail/locus_arrangement.hpp>
#include "math/fractional.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <iterator>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace cppcrystal::group {

namespace {

// Number of free coordinates of a Wyckoff position = rank of its representative
// coordinate operator (the projector onto the canonical coordinate).
[[nodiscard]] int rank_of(Matrix3i const &rotation) {
  Eigen::FullPivLU<Matrix3d> lu(rotation.cast<double>());
  lu.setThreshold(0.5); // integer entries: cleanly separates rank levels
  return static_cast<int>(lu.rank());
}

} // namespace

WyckoffPosition SpaceGroup::build_position(data::WyckoffEntry const &entry,
                                           Operations const &conv_ops) {
  data::WyckoffCoordinate const wc =
      data::wyckoff_coordinate(entry.global_index);

  // A generic seed, projected onto the position, lands at a generic point whose
  // stabilizer is exactly the site-symmetry group (for a 0-DOF position the
  // projector ignores the seed and pins the fixed point). Exact (rational)
  // database coordinates are compared at the default symprec.
  Vector3d const seed{0.1357, 0.2468, 0.3791};
  Vector3d const p0 = math::wrap_to_unit_cell(
      wc.rotation.cast<double>() * seed + wc.translation);

  auto partition = detail::partition_orbit(
      conv_ops, p0, [](Vector3d const &p) { return math::wrap_to_unit_cell(p); },
      [](Vector3d const &a, Vector3d const &b) {
        return math::same_point(a, b, kDefaultSymprec);
      });

  return WyckoffPosition{wc.multiplicity,
                         entry.letter,
                         rank_of(wc.rotation),
                         data::site_symmetry_symbol(entry.global_index),
                         wc.rotation,
                         wc.translation,
                         std::move(partition.orbit_ops),
                         std::move(partition.site_symmetry)};
}

Result<SpaceGroup> SpaceGroup::from_hall_number(int hall_number) {
  data::SpacegroupType const type = data::spacegroup_type(hall_number);
  if (type.number == 0) {
    return leaf::new_error(
        e_message{"SpaceGroup::from_hall_number: Hall number out of range "
                  "(expected 1..530)"});
  }

  SpaceGroup sg;
  sg.type_ = type;
  sg.operations_ = data::operations_from_database(hall_number);

  std::vector<data::WyckoffEntry> const entries =
      data::wyckoff_entries(hall_number);
  sg.positions_.reserve(entries.size());
  std::ranges::transform(entries, std::back_inserter(sg.positions_),
                         [&](data::WyckoffEntry const &entry) {
                           return build_position(entry, sg.operations_);
                         });
  return sg;
}

Result<SpaceGroup> SpaceGroup::from_layer_hall(int hall_number) {
  data::SpacegroupType const type = data::spacegroup_type(hall_number);
  if (hall_number >= 0 || type.number == 0) {
    return leaf::new_error(e_message{
        "SpaceGroup::from_layer_hall: layer hall number out of range "
        "(expected -1..-116)"});
  }

  SpaceGroup lg;
  lg.type_ = type;
  lg.operations_ = data::operations_from_database(hall_number);

  std::vector<data::WyckoffEntry> const entries =
      data::wyckoff_entries(hall_number);
  lg.positions_.reserve(entries.size());
  std::ranges::transform(
      entries, std::back_inserter(lg.positions_),
      [&](const auto &entry) { return build_position(entry, lg.operations_); });
  return lg;
}

Result<SpaceGroup> SpaceGroup::from_layer_number(int layer_number) {
  // The layer-group number -> default-Hall map is resolved at compile time.
  int const hall = data::layer_default_hall(layer_number);
  if (hall == 0) {
    return leaf::new_error(e_message{
        "SpaceGroup::from_layer_number: layer-group number out of range "
        "(expected 1..80)"});
  }
  return from_layer_hall(hall);
}

Result<SpaceGroup> SpaceGroup::from_number(int spacegroup_number) {
  std::span<int const> const halls =
      data::halls_with_number(spacegroup_number);
  if (halls.empty()) {
    return leaf::new_error(
        e_message{"SpaceGroup::from_number: international number out of range "
                  "(expected 1..230)"});
  }
  return from_hall_number(halls.front());
}

Result<WyckoffPosition const *> SpaceGroup::wyckoff(char letter) const {
  if (auto const it = std::ranges::find(positions_, letter,
                                        &WyckoffPosition::letter);
      it != positions_.end()) {
    return &*it;
  }
  return leaf::new_error(
      e_message{std::string("SpaceGroup::wyckoff: no Wyckoff position '") +
                letter + "' in this space group"});
}

} // namespace cppcrystal::group
