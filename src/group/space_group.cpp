#include <cppcrystal/group/space_group.hpp>

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/group/detail/locus_arrangement.hpp>
#include "math/fractional.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <mutex>
#include <shared_mutex>
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

SpaceGroup::SpaceGroup(HallNumber hall)
    : hall_(hall), operations_(data::operations_from_database(hall)) {
  std::vector<data::WyckoffEntry> const entries = data::wyckoff_entries(hall);
  positions_.reserve(entries.size());
  std::ranges::transform(entries, std::back_inserter(positions_),
                         [&](data::WyckoffEntry const &entry) {
                           return build_position(entry, operations_);
                         });
}

// The Flyweight store: one slot per setting per family, filled on first use.
// A shared_mutex keeps the common case (an already-built setting) to a shared
// lock, which matters because SpaceGroup::of sits under every generation call.
SpaceGroup const &SpaceGroup::of(HallNumber hall) {
  struct Store {
    std::shared_mutex mutex;
    std::array<std::unique_ptr<SpaceGroup>, kSpaceHallSettings> space;
    std::array<std::unique_ptr<SpaceGroup>, kLayerHallSettings> layer;

    [[nodiscard]] std::unique_ptr<SpaceGroup> &slot(HallNumber h) noexcept {
      auto const i = static_cast<std::size_t>(h.index()) - 1;
      return h.family() == GroupFamily::layer ? layer[i] : space[i];
    }
  };
  static Store store;

  {
    std::shared_lock const read(store.mutex);
    if (auto const &built = store.slot(hall)) {
      return *built;
    }
  }
  // Build outside the lock: construction is independent per setting, and a
  // duplicate build is cheaper than serialising every miss.
  auto fresh = std::unique_ptr<SpaceGroup>(new SpaceGroup(hall));
  std::unique_lock const write(store.mutex);
  auto &slot = store.slot(hall);
  if (!slot) {
    slot = std::move(fresh);
  }
  return *slot;
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
