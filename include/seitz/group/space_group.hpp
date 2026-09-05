#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/group/wyckoff.hpp>

#include <span>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::group {

// A space group as a standalone, structure-free object: from a Hall (or
// international) number, owning its conventional operations and its Wyckoff
// positions. No atomic structure involved -- the queryable face of the symmetry
// database, under crystal generation and group-relation work.
class SpaceGroup : public GroupBase {
public:
  // The group of a Hall setting, as a Boost.Flyweight: one immutable object
  // per setting, built on first use, so Wyckoff construction is paid once per
  // setting, not per query. Total -- a HallNumber always names a setting that
  // exists. Layer groups come through the same door, family carried by the
  // key. Race-free; warmup() primes it.
  [[nodiscard]] static SpaceGroup const &of(HallNumber hall);

  // A private, unshared instance; the flyweight factory builds settings
  // through it, callers want `of`.
  explicit SpaceGroup(HallNumber hall);

  // The group of an international number, in its default (first) Hall setting.
  // Errors if the number is out of range for the family (1..230 for space
  // groups, 1..80 for layer groups).
  [[nodiscard]] static Result<SpaceGroup const *>
  from_number(GroupFamily family, int number);

  [[nodiscard]] HallNumber hall() const noexcept { return hall_; }
  [[nodiscard]] data::SpacegroupType const &type() const noexcept {
    return data::spacegroup_type(hall_);
  }
  [[nodiscard]] data::Centering centering() const noexcept {
    return type().centering;
  }

private:
  // One Wyckoff position, by partitioning the conventional operations into its
  // orbit (coset representatives) and its stabilizer. A member so it can reach
  // Wyckoff's private constructor.
  [[nodiscard]] static Wyckoff build_position(int global_index, int letter,
                                              Operations const &conv_ops);

  HallNumber hall_;
};

} // namespace seitz::group

#pragma GCC visibility pop
