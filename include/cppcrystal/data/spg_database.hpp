#pragma once

#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/catalog.hpp>
#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/spacegroup_metadata_tables.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <string_view>

// Access to the built-in space-group database. The raw data tables are
// generated into two headers: spacegroup_metadata_tables.hpp (plain metadata,
// included here to back the constexpr catalog) and
// spacegroup_operation_tables.hpp (the encoded operations — a
// compile-time-only compaction, included only by the .cpp). The metadata is
// decoded once at compile time into a Catalog per family; the operations go
// through a lazily-built cache, since they carry Eigen, which is not a literal
// type.
// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal::data {

// Bravais centering.
enum class Centering {
  error,
  primitive,
  body,
  face,
  a_face,
  b_face,
  c_face,
  base,
  r_center,
};

// One Hall setting's metadata. The Hall number itself is NOT a field: it is
// the key that addresses the row, and storing it alongside would be a second
// source of truth.
struct SpacegroupType {
  int number = 0; // international space-group / layer-group number
  std::string_view schoenflies;
  std::string_view hall_symbol;
  std::string_view international;
  std::string_view international_full;
  std::string_view international_short;
  std::string_view choice;
  Centering centering = Centering::error;
  int pointgroup_number = 0; // 1..32
};

inline constexpr int kNumSpacegroups = 230;
inline constexpr int kNumLayerGroups = 80;
inline constexpr int kNumPointgroups = 32;

// ---- family policies -------------------------------------------------------

// The two space-group families differ only in which generated table they
// decode and how many settings and groups they have.
template <GroupFamily F> struct SpacegroupFamily {
  using Row = SpacegroupType;
  using Key = HallNumber;

  static constexpr std::size_t count =
      static_cast<std::size_t>(hall_settings(F));
  static constexpr int groups =
      F == GroupFamily::layer ? kNumLayerGroups : kNumSpacegroups;

  [[nodiscard]] static constexpr int index_of(HallNumber key) noexcept {
    return key.index();
  }

  [[nodiscard]] static constexpr auto const &raw() noexcept {
    if constexpr (F == GroupFamily::layer) {
      return kLayerGroupTypes;
    } else {
      return kSpacegroupTypes;
    }
  }

  // The generated tables are 1-based with a dummy row 0; the catalog is not.
  [[nodiscard]] static constexpr std::array<Row, count> decode() {
    std::array<Row, count> rows{};
    for (std::size_t i = 0; i < count; ++i) {
      auto const &r = raw()[i + 1];
      rows[i] = SpacegroupType{r.number,
                               r.schoenflies,
                               r.hall_symbol,
                               r.international,
                               r.international_full,
                               r.international_short,
                               r.choice,
                               static_cast<Centering>(r.centering),
                               r.pointgroup_number};
    }
    return rows;
  }
};

using SpaceGroupFamily = SpacegroupFamily<GroupFamily::space>;
using LayerGroupFamily = SpacegroupFamily<GroupFamily::layer>;

// ---- keyed access ----------------------------------------------------------

// The setting a Hall number names. Total: HallNumber has no invalid state.
template <GroupFamily F>
[[nodiscard]] constexpr SpacegroupType const &
spacegroup_type(HallNumber hall) noexcept {
  return kCatalog<SpacegroupFamily<F>>[hall];
}

[[nodiscard]] constexpr SpacegroupType const &
spacegroup_type(HallNumber hall) noexcept {
  return hall.family() == GroupFamily::layer
             ? spacegroup_type<GroupFamily::layer>(hall)
             : spacegroup_type<GroupFamily::space>(hall);
}

// ---- inverted indices ------------------------------------------------------

// international number -> its Hall settings, ascending (so the first is the
// default setting). Built at compile time from the catalog.
template <GroupFamily F>
inline constexpr auto kHallsByNumber =
    detail::bucket_index<static_cast<std::size_t>(hall_settings(F)),
                         static_cast<std::size_t>(SpacegroupFamily<F>::groups) +
                             1>([](int index) {
      return kCatalog<SpacegroupFamily<F>>.rows[
              static_cast<std::size_t>(index) - 1].number;
    });

// Every Hall setting index of an international number; empty if out of range.
template <GroupFamily F>
[[nodiscard]] constexpr std::span<int const>
halls_with_number(int number) noexcept {
  return kHallsByNumber<F>[number];
}

// The default (first) Hall setting of an international number.
template <GroupFamily F>
[[nodiscard]] constexpr std::optional<HallNumber>
default_hall(int number) noexcept {
  auto const halls = halls_with_number<F>(number);
  return halls.empty() ? std::nullopt : HallNumber::of(F, halls.front());
}

// point-group number -> the default Hall setting index of every group with
// that point group, ascending by group number. The candidate list of the
// space-group search: a found point group narrows the candidates to these.
template <GroupFamily F>
inline constexpr auto kDefaultHallsByPointgroup =
    detail::bucket_index<static_cast<std::size_t>(SpacegroupFamily<F>::groups),
                         static_cast<std::size_t>(kNumPointgroups) + 1>(
        [](int number) {
          return kCatalog<SpacegroupFamily<F>>
              .rows[static_cast<std::size_t>(
                        kHallsByNumber<F>[number].front()) - 1]
              .pointgroup_number;
        },
        [](int number) { return kHallsByNumber<F>[number].front(); });

template <GroupFamily F>
[[nodiscard]] constexpr std::span<int const>
default_halls_with_pointgroup(int pointgroup_number) noexcept {
  return kDefaultHallsByPointgroup<F>[pointgroup_number];
}

// ---- compile-time integrity guards -----------------------------------------

namespace detail {
// Every group lands in exactly one point-group bucket, under its own number.
template <GroupFamily F> constexpr bool pointgroup_index_well_formed() {
  std::size_t total = 0;
  for (int pg = 0; pg <= kNumPointgroups; ++pg) {
    for (int index : kDefaultHallsByPointgroup<F>[pg]) {
      if (kCatalog<SpacegroupFamily<F>>
              .rows[static_cast<std::size_t>(index) - 1]
              .pointgroup_number != pg) {
        return false;
      }
      ++total;
    }
  }
  return total == static_cast<std::size_t>(SpacegroupFamily<F>::groups);
}

// The settings must cover the group numbers in ascending blocks, so the first
// setting of each number is its canonical one. Breaks the build if the
// generated tables ever drift.
template <GroupFamily F> constexpr bool catalog_well_formed() {
  auto const &rows = kCatalog<SpacegroupFamily<F>>.rows;
  constexpr int groups = SpacegroupFamily<F>::groups;
  return SpacegroupFamily<F>::raw().size() == SpacegroupFamily<F>::count + 1 &&
         std::ranges::all_of(rows,
                             [](SpacegroupType const &t) {
                               return t.number >= 1 && t.number <= groups;
                             }) &&
         std::ranges::is_sorted(rows, {}, &SpacegroupType::number) &&
         std::ranges::max(rows, {}, &SpacegroupType::number).number == groups &&
         std::ranges::all_of(std::views::iota(1, groups + 1), [](int number) {
           return default_hall<F>(number).has_value();
         });
}
} // namespace detail

static_assert(detail::pointgroup_index_well_formed<GroupFamily::space>());
static_assert(detail::pointgroup_index_well_formed<GroupFamily::layer>());
static_assert(detail::catalog_well_formed<GroupFamily::space>());
static_assert(detail::catalog_well_formed<GroupFamily::layer>());

static_assert(default_halls_with_pointgroup<GroupFamily::space>(0).empty());
// m-3m: 221..230
static_assert(default_halls_with_pointgroup<GroupFamily::space>(32).size() ==
              10);
// P2: unique axis b, c, a
static_assert(halls_with_number<GroupFamily::space>(3).size() == 3);
static_assert(default_hall<GroupFamily::space>(1)->index() == 1);
static_assert(!default_hall<GroupFamily::space>(231).has_value());
// p1 is the first layer setting
static_assert(default_hall<GroupFamily::layer>(1)->index() == 1);
static_assert(default_hall<GroupFamily::layer>(80).has_value());
static_assert(!default_hall<GroupFamily::layer>(81).has_value());
static_assert(!HallNumber::of(GroupFamily::space, 0).has_value());
static_assert(
    !HallNumber::of(GroupFamily::space, kSpaceHallSettings + 1).has_value());
static_assert(
    !HallNumber::of(GroupFamily::layer, kLayerHallSettings + 1).has_value());
static_assert(spacegroup_type(*HallNumber::of(GroupFamily::layer, 1)).number ==
              1);
static_assert(spacegroup_type(*HallNumber::of(GroupFamily::layer,
                                              kLayerHallSettings))
                  .number == kNumLayerGroups);

// ---- operations ------------------------------------------------------------

// The symmetry operations of a Hall setting, materialised once per setting.
[[nodiscard]] Operations const &operations_from_database(HallNumber hall);

} // namespace cppcrystal::data

#pragma GCC visibility pop
