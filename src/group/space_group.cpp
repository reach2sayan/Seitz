#include <spglib/group/space_group.hpp>

#include <spglib/data/sitesym_database.hpp>
#include <spglib/math/fractional.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace spglib::group {

namespace {
constexpr double kGroupTol = 1e-5;

[[nodiscard]] bool same_point(Vector3d const &a, Vector3d const &b) {
  return math::frac_displacement(a, b).cwiseAbs().maxCoeff() < kGroupTol;
}

// Number of free coordinates of a Wyckoff position = rank of its representative
// coordinate operator (the projector onto the canonical coordinate).
[[nodiscard]] int rank_of(Matrix3i const &rotation) {
  Eigen::FullPivLU<Matrix3d> lu(rotation.cast<double>());
  lu.setThreshold(0.5); // integer entries: cleanly separates rank levels
  return static_cast<int>(lu.rank());
}

} // namespace

WyckoffPosition SpaceGroup::build_position(data::WyckoffEntry const &entry,
                                           SymmetryOperations const &conv_ops) {
  data::WyckoffCoordinate const wc =
      data::wyckoff_coordinate(entry.global_index);

  // A generic seed, projected onto the position, lands at a generic point whose
  // stabilizer is exactly the site-symmetry group (for a 0-DOF position the
  // projector ignores the seed and pins the fixed point).
  Vector3d const seed{0.1357, 0.2468, 0.3791};
  Vector3d const p0 = math::wrap_to_unit_cell(
      wc.rotation.cast<double>() * seed + wc.translation);

  SymmetryOperations orbit_ops; // coset representatives, one per orbit point
  SymmetryOperations site_symmetry; // stabilizer of p0
  std::vector<Vector3d> distinct_images;
  for (auto const &op : conv_ops) {
    Vector3d const image = math::wrap_to_unit_cell(op.apply(p0));
    if (same_point(image, p0)) {
      site_symmetry.push_back(op);
    }
    bool const seen =
        std::ranges::any_of(distinct_images, [&](Vector3d const &q) {
          return same_point(q, image);
        });
    if (!seen) {
      distinct_images.push_back(image);
      orbit_ops.push_back(op);
    }
  }

  return WyckoffPosition{
      wc.multiplicity,      entry.letter,
      rank_of(wc.rotation), data::site_symmetry_symbol(entry.global_index),
      wc.rotation,          wc.translation,
      std::move(orbit_ops), std::move(site_symmetry)};
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
  for (auto const &entry : entries) {
    sg.positions_.push_back(build_position(entry, sg.operations_));
  }
  return sg;
}

Result<SpaceGroup> SpaceGroup::from_number(int spacegroup_number) {
  std::span<int const> const halls =
      data::kCatalog.with_number(spacegroup_number);
  if (halls.empty()) {
    return leaf::new_error(
        e_message{"SpaceGroup::from_number: international number out of range "
                  "(expected 1..230)"});
  }
  return from_hall_number(halls.front());
}

Result<WyckoffPosition const *> SpaceGroup::wyckoff(char letter) const {
  for (auto const &wp : positions_) {
    if (wp.letter() == letter) {
      return &wp;
    }
  }
  return leaf::new_error(
      e_message{std::string("SpaceGroup::wyckoff: no Wyckoff position '") +
                letter + "' in this space group"});
}

} // namespace spglib::group
