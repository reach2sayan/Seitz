#pragma once

#include <compare>
#include <optional>

#pragma GCC visibility push(default)

namespace seitz {

// The two crystallographic families the determination pipeline handles. A
// layer group is 2D-periodic with its own 116-setting Hall database; the rest
// of the search is shared with the 3D space groups, so the family is carried as
// a compile-time parameter and every family-dependent branch is `if constexpr`.
enum class GroupFamily { space, layer };

// Number of Hall settings per family. Static-asserted against the generated
// tables in data/spg_database.hpp.
inline constexpr int kSpaceHallSettings = 530;
inline constexpr int kLayerHallSettings = 116;

[[nodiscard]] constexpr int hall_settings(GroupFamily family) noexcept {
  return family == GroupFamily::layer ? kLayerHallSettings : kSpaceHallSettings;
}

// A validated Hall-setting key: a family plus a 1-based index into that
// family's settings. Replaces the old signed convention (a 3D setting as
// 1..530, a layer setting as -1..-116), so a catalog lookup keyed by one of
// these is total and the catalogs need no sentinel row.
//
// There is no default constructor and no invalid state: `of` is the only way
// in, and it rejects an out-of-range index.
class HallNumber {
public:
  [[nodiscard]] static constexpr std::optional<HallNumber>
  of(GroupFamily family, int index) noexcept {
    if (index < 1 || index > hall_settings(family)) {
      return std::nullopt;
    }
    return HallNumber{family, index};
  }

  [[nodiscard]] constexpr GroupFamily family() const noexcept {
    return family_;
  }
  [[nodiscard]] constexpr int index() const noexcept { return index_; }

  [[nodiscard]] friend constexpr bool operator==(HallNumber,
                                                 HallNumber) = default;
  [[nodiscard]] friend constexpr std::strong_ordering
  operator<=>(HallNumber, HallNumber) = default;

private:
  constexpr HallNumber(GroupFamily family, int index) noexcept
      : family_{family}, index_{index} {}

  GroupFamily family_;
  int index_;
};

// Number of magnetic space-group (UNI) numbers.
inline constexpr int kUniNumbers = 1651;

// A validated UNI (magnetic space-group) number, 1..1651.
class UniNumber {
public:
  [[nodiscard]] static constexpr std::optional<UniNumber>
  of(int number) noexcept {
    if (number < 1 || number > kUniNumbers) {
      return std::nullopt;
    }
    return UniNumber{number};
  }

  [[nodiscard]] constexpr int value() const noexcept { return value_; }

  [[nodiscard]] friend constexpr bool operator==(UniNumber,
                                                 UniNumber) = default;
  [[nodiscard]] friend constexpr std::strong_ordering
  operator<=>(UniNumber, UniNumber) = default;

private:
  explicit constexpr UniNumber(int value) noexcept : value_{value} {}
  int value_;
};

} // namespace seitz

#pragma GCC visibility pop
