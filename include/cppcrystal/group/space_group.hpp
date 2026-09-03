#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/group/wyckoff.hpp>

#include <span>
#include <vector>

namespace cppcrystal::group {

// A space group as a standalone, structure-free object: built from a Hall
// number (or international number), it owns its conventional symmetry
// operations and its Wyckoff positions as first-class objects. No atomic
// structure is involved — this is the queryable catalog face of the symmetry
// database, the foundation for crystal generation and group-relation work.
//
// Wyckoff positions are exposed as a span of objects rather than parallel
// arrays of letters / multiplicities / symbols.
class SpaceGroup : public GroupBase {
public:
  // The group of a Hall setting. A Flyweight: one immutable object per setting,
  // built on first use and shared thereafter, so the Wyckoff construction is
  // paid once per setting rather than once per query. Total — a HallNumber
  // cannot name a setting that does not exist. Layer groups come through the
  // same door, with the family carried by the key.
  //
  // Thread-safety: the per-setting cache is guarded, so concurrent first-calls
  // are safe. warmup() primes it.
  [[nodiscard]] static SpaceGroup const &of(HallNumber hall);

  // The group of an international number, in its default (first) Hall setting.
  // Errors if the number is out of range for the family (1..230 for space
  // groups, 1..80 for layer groups).
  [[nodiscard]] static Result<SpaceGroup const *> from_number(GroupFamily family,
                                                              int number);

  [[nodiscard]] HallNumber hall() const noexcept { return hall_; }
  [[nodiscard]] data::SpacegroupType const &type() const noexcept {
    return data::spacegroup_type(hall_);
  }
  [[nodiscard]] data::Centering centering() const noexcept {
    return type().centering;
  }

private:
  explicit SpaceGroup(HallNumber hall);

  // Build one Wyckoff position by partitioning the conventional operations into
  // its orbit (coset representatives) and site-symmetry stabilizer. A member so
  // it can reach Wyckoff's private constructor (SpaceGroup is a
  // friend).
  [[nodiscard]] static Wyckoff
  build_position(data::WyckoffEntry const &entry,
                 Operations const &conv_ops);

  HallNumber hall_;
};

} // namespace cppcrystal::group
