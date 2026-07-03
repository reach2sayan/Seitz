#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/group/wyckoff_position.hpp>

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
class SpaceGroup {
public:
  // Named factories.
  [[nodiscard]] static Result<SpaceGroup> from_hall_number(int hall_number);
  // Build from an international number (1..230), choosing the default (first)
  // Hall setting for that number.
  [[nodiscard]] static Result<SpaceGroup> from_number(int spacegroup_number);

  // Build a *layer group* (2D-periodic) as a structure-free object. Layer
  // groups use a negative-Hall-number convention (settings -1..-116 for the 80
  // layer groups), so `number()` then returns the layer-group number (1..80)
  // and `hall_number()` the negative setting. `from_layer_hall` takes a
  // negative hall (-1..-116); `from_layer_number` takes a layer-group number
  // (1..80) and picks its first setting.
  [[nodiscard]] static Result<SpaceGroup> from_layer_hall(int hall_number);
  [[nodiscard]] static Result<SpaceGroup> from_layer_number(int layer_number);

  [[nodiscard]] int hall_number() const noexcept { return type_.hall_number; }
  [[nodiscard]] int number() const noexcept { return type_.number; }
  [[nodiscard]] data::SpacegroupType const &type() const noexcept {
    return type_;
  }
  [[nodiscard]] std::string_view international_symbol() const noexcept {
    return type_.international_short;
  }
  [[nodiscard]] data::Centering centering() const noexcept {
    return type_.centering;
  }

  // The conventional symmetry operations of the setting.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return operations_;
  }

  // The Wyckoff positions, ordered by ascending letter (a, b, c, ...). Owned as
  // objects, not parallel arrays.
  [[nodiscard]] std::span<WyckoffPosition const> wyckoffs() const noexcept {
    return positions_;
  }

  // Look up a Wyckoff position by letter ('a' = the most special). Errors if
  // the letter is out of range for this group.
  [[nodiscard]] Result<WyckoffPosition const *> wyckoff(char letter) const;

private:
  SpaceGroup() = default;

  // Build one Wyckoff position by partitioning the conventional operations into
  // its orbit (coset representatives) and site-symmetry stabilizer. A member so
  // it can reach WyckoffPosition's private constructor (SpaceGroup is a
  // friend).
  [[nodiscard]] static WyckoffPosition
  build_position(data::WyckoffEntry const &entry,
                 SymmetryOperations const &conv_ops);

  data::SpacegroupType type_;
  SymmetryOperations operations_;
  std::vector<WyckoffPosition> positions_;
};

} // namespace cppcrystal::group
