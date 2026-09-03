#include <cppcrystal/group/space_group.hpp>

#include "data/sitesym_database.hpp"
#include "group/locus_arrangement.hpp"
#include "math/fractional.hpp"
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>

#include <Eigen/Dense>

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <boost/flyweight/set_factory.hpp>

#include <algorithm>
#include <iterator>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace cppcrystal::group {

namespace {

// Wyckoff letters run a..z then A.. for the rare groups with more than 26.
[[nodiscard]] char letter_of(int index) {
  return index < 26 ? static_cast<char>('a' + index)
                    : static_cast<char>('A' + (index - 26));
}

// The locus of a Wyckoff position, read off its representative coordinate
// operator (which IS the projector onto the position): the number of free
// coordinates is the operator's rank and its image spans the locus directions.
// Both come from one decomposition, so the basis has exactly `dof` columns --
// deriving them separately (an LU rank against an SVD column space) can
// disagree on a near-threshold singular value and silently mis-size the basis.
struct Locus {
  int dof = 0;
  Matrix3d basis = Matrix3d::Zero(); // columns 0..dof-1
};

[[nodiscard]] Locus locus_of(Matrix3i const &rotation) {
  Matrix3d const projector = rotation.cast<double>();
  Eigen::FullPivLU<Matrix3d> lu(projector);
  lu.setThreshold(0.5); // integer entries: cleanly separates rank levels
  Locus out;
  out.dof = static_cast<int>(lu.rank());
  if (out.dof > 0) {
    out.basis.leftCols(out.dof) = lu.image(projector);
  }
  return out;
}

} // namespace

Wyckoff SpaceGroup::build_position(int global_index, int letter,
                                   Operations const &conv_ops) {
  data::WyckoffCoordinate const wc = data::wyckoff_coordinate(global_index);

  // A generic seed, projected onto the position, lands at a generic point whose
  // stabilizer is exactly the site-symmetry group (for a 0-DOF position the
  // projector ignores the seed and pins the fixed point). Exact (rational)
  // database coordinates are compared at the default symprec.
  Vector3d const seed{0.1357, 0.2468, 0.3791};
  Vector3d const p0 = math::wrap_to_unit_cell(
      wc.rotation.cast<double>() * seed + wc.translation);

  auto partition = detail::partition_orbit(
      conv_ops, p0,
      [](Vector3d const &p) { return math::wrap_to_unit_cell(p); },
      [](Vector3d const &a, Vector3d const &b) {
        return math::same_point(a, b, kDefaultSymprec);
      });

  Locus const locus = locus_of(wc.rotation);
  return Wyckoff{wc.multiplicity,
                 locus.dof,
                 letter_of(letter),
                 data::site_symmetry_symbol(global_index),
                 wc.translation,
                 locus.basis,
                 wc.rotation.cast<double>(),
                 wc.translation,
                 std::move(partition.orbit_ops),
                 std::move(partition.site_symmetry)};
}

SpaceGroup::SpaceGroup(HallNumber hall) : hall_(hall) {
  data::SpacegroupType const &t = data::spacegroup_type(hall);
  number_ = t.number;
  symbol_ = t.international_short;
  operations_ = data::operations_from_database(hall);

  std::vector<data::WyckoffEntry> const entries = data::wyckoff_entries(hall);
  positions_.reserve(entries.size());
  std::ranges::transform(entries, std::back_inserter(positions_),
                         [&](data::WyckoffEntry const &entry) {
                           return build_position(entry.global_index,
                                                 entry.letter, operations_);
                         });
}

// The Flyweight (Boost.Flyweight): one immutable SpaceGroup per HallNumber,
// built from the key on first use and shared thereafter. `set_factory` keys on
// HallNumber's ordering (no hash needed) and keeps addresses stable;
// `no_tracking` pins every built setting for the program's lifetime, which is
// what lets `of` hand out plain references; the default locking policy makes
// concurrent first-calls safe. warmup() primes it.
SpaceGroup const &SpaceGroup::of(HallNumber hall) {
  using Shared =
      boost::flyweight<boost::flyweights::key_value<HallNumber, SpaceGroup>,
                       boost::flyweights::set_factory<>,
                       boost::flyweights::no_tracking>;
  return Shared(hall).get();
}

Result<SpaceGroup const *> SpaceGroup::from_number(GroupFamily family,
                                                   int number) {
  auto const hall = family == GroupFamily::layer
                        ? data::default_hall<GroupFamily::layer>(number)
                        : data::default_hall<GroupFamily::space>(number);
  if (!hall) {
    return leaf::new_error(e_message{
        family == GroupFamily::layer
            ? "SpaceGroup::from_number: layer-group number out of range "
              "(expected 1..80)"
            : "SpaceGroup::from_number: international number out of range "
              "(expected 1..230)"});
  }
  return &of(*hall);
}

} // namespace cppcrystal::group
