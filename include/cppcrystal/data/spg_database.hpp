#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/spacegroup_metadata_tables.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <ranges>
#include <span>
#include <string_view>

// Access to the built-in space-group database (530 Hall settings). The raw data
// tables are generated into two headers: spacegroup_metadata_tables.hpp (plain
// metadata, included here to back the constexpr catalog) and
// spacegroup_operation_tables.hpp (the encoded operations — a
// compile-time-only compaction, included only by the .cpp). Here we decode the
// metadata once at compile time into a constexpr multi-keyed catalog, and expose
// the symmetry operations through a lazily-built cache (the ops carry Eigen,
// which is not a literal type, so they cannot live in a constexpr global). Layer
// groups (negative Hall numbers -1..-116) share this machinery and are decoded
// into the parallel constexpr kLayerCatalog below.
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

struct SpacegroupType {
  int hall_number = 0; // 1..530 (the catalog's unique key)
  int number = 0;      // international space-group number 1..230
  std::string_view schoenflies;
  std::string_view hall_symbol;
  std::string_view international;
  std::string_view international_full;
  std::string_view international_short;
  std::string_view choice;
  Centering centering = Centering::error;
  int pointgroup_number = 0; // 1..32
};

inline constexpr int kNumHallNumbers = 530;

// A compile-time, multi-keyed catalog of the 530 Hall settings. Lookup by Hall
// number is a direct index (Hall numbers are the contiguous range 1..530);
// lookup by international number uses a sorted index permutation with an
// equal_range-style query. Index 0 is a zero-valued sentinel returned for
// out-of-range Hall numbers.
struct SpacegroupCatalog {
  std::array<SpacegroupType, kNumHallNumbers + 1> by_hall{};
  // Hall numbers 1..530 permuted so that by_hall[halls_by_number[i]].number is
  // non-decreasing (ties ordered by Hall number).
  std::array<int, kNumHallNumbers> halls_by_number{};

  // The setting for a Hall number (1..530); the sentinel (number 0) otherwise.
  [[nodiscard]] constexpr SpacegroupType const &at(int hall) const noexcept {
    return detail::at_or_sentinel(by_hall, hall);
  }

  // Every Hall setting of a given international number, as a span of Hall
  // numbers (empty if none).
  [[nodiscard]] constexpr std::span<int const>
  with_number(int number) const noexcept {
    auto const less_n = [this](int hall, int n) {
      return by_hall[static_cast<std::size_t>(hall)].number < n;
    };
    auto const n_less = [this](int n, int hall) {
      return n < by_hall[static_cast<std::size_t>(hall)].number;
    };
    auto const first = std::lower_bound(halls_by_number.begin(),
                                        halls_by_number.end(), number, less_n);
    auto const last = std::upper_bound(halls_by_number.begin(),
                                       halls_by_number.end(), number, n_less);
    return {halls_by_number.data() + (first - halls_by_number.begin()),
            static_cast<std::size_t>(last - first)};
  }
};

// The catalog, decoded from the raw tables at compile time.
inline constexpr SpacegroupCatalog kCatalog = [] {
  SpacegroupCatalog c{};
  for (std::size_t h = 0; h < kSpacegroupTypes.size(); ++h) {
    auto const &r = kSpacegroupTypes[h];
    c.by_hall[h] = SpacegroupType{static_cast<int>(h),
                                  r.number,
                                  r.schoenflies,
                                  r.hall_symbol,
                                  r.international,
                                  r.international_full,
                                  r.international_short,
                                  r.choice,
                                  static_cast<Centering>(r.centering),
                                  r.pointgroup_number};
  }
  for (int i = 0; i < kNumHallNumbers; ++i) {
    c.halls_by_number[static_cast<std::size_t>(i)] = i + 1;
  }
  std::sort(c.halls_by_number.begin(), c.halls_by_number.end(),
            [&c](int a, int b) {
              auto const &ta = c.by_hall[static_cast<std::size_t>(a)];
              auto const &tb = c.by_hall[static_cast<std::size_t>(b)];
              return ta.number != tb.number ? ta.number < tb.number : a < b;
            });
  return c;
}();

// Number of layer-group Hall settings (80 layer groups, 116 settings). Layer
// groups use a negative-hall-number convention (settings -1..-116), so they
// share the symmetry-operation array and the bulk of the matching machinery with
// the 3D space groups.
inline constexpr int kNumLayerHallNumbers = 116;

// Number of distinct layer groups (the 116 settings collapse onto 80 numbers).
inline constexpr int kNumLayerGroups = 80;

[[nodiscard]] constexpr int num_layer_groups() noexcept {
  return kNumLayerGroups;
}

// Whether `hall_number` names a layer-group Hall setting (-1..-116), the layer
// analogue of a positive 3D Hall number being in 1..530.
[[nodiscard]] constexpr bool layer_hall_in_range(int hall_number) noexcept {
  return hall_number <= -1 && hall_number >= -kNumLayerHallNumbers;
}

// Whether `layer_number` is a valid layer-group number (1..80).
[[nodiscard]] constexpr bool layer_number_in_range(int layer_number) noexcept {
  return layer_number >= 1 && layer_number <= kNumLayerGroups;
}

// The layer-group analogue of SpacegroupCatalog, indexed by the negation of the
// (negative) layer hall number: by_neg_hall[-hall]. Index 0 is the sentinel.
struct LayerCatalog {
  std::array<SpacegroupType, kNumLayerHallNumbers + 1> by_neg_hall{};

  // The setting for a layer hall number (-1..-116); the sentinel otherwise.
  [[nodiscard]] constexpr SpacegroupType const &at(int hall) const noexcept {
    return detail::at_or_sentinel(by_neg_hall, -hall);
  }
};

inline constexpr LayerCatalog kLayerCatalog = [] {
  LayerCatalog c{};
  for (std::size_t i = 0; i < kLayerGroupTypes.size(); ++i) {
    auto const &r = kLayerGroupTypes[i];
    c.by_neg_hall[i] = SpacegroupType{-static_cast<int>(i),
                                      r.number,
                                      r.schoenflies,
                                      r.hall_symbol,
                                      r.international,
                                      r.international_full,
                                      r.international_short,
                                      r.choice,
                                      static_cast<Centering>(r.centering),
                                      r.pointgroup_number};
  }
  return c;
}();

// The setting for a Hall number; a zero-valued SpacegroupType (number 0) if out
// of range. Positive hall numbers (1..530) select 3D space groups; negative
// hall numbers (-1..-116) select layer-group settings. constexpr — usable in
// constant expressions.
[[nodiscard]] constexpr SpacegroupType
spacegroup_type(int hall_number) noexcept {
  return detail::hall_indexed(kCatalog.by_hall, kLayerCatalog.by_neg_hall,
                              hall_number);
}

// number -> first (default) layer Hall setting, built once at compile time.
// Walking the settings from the last down leaves each number's smallest
// (first/canonical) setting in place; entry 0 stays 0 for out-of-range.
inline constexpr auto kLayerDefaultHall = [] {
  std::array<int, kNumLayerGroups + 1> t{};
  for (int neg = kNumLayerHallNumbers; neg >= 1; --neg) {
    int const number = kLayerCatalog.by_neg_hall[static_cast<std::size_t>(neg)].number;
    if (number >= 1 && number <= kNumLayerGroups) {
      t[static_cast<std::size_t>(number)] = -neg;
    }
  }
  return t;
}();

// The first (default) layer Hall setting (-1..-116) for a layer-group number
// (1..80); 0 if out of range.
[[nodiscard]] constexpr int layer_default_hall(int layer_number) noexcept {
  return layer_number_in_range(layer_number)
             ? kLayerDefaultHall[static_cast<std::size_t>(layer_number)]
             : 0;
}

// Compile-time integrity guards on the layer catalog: the settings must cover
// the layer-group numbers 1..80 in ascending blocks, with every number resolving
// to a default Hall setting. These break the build if the generated layer tables
// ever drift (the layer counterpart of hall_classification.hpp's guards).
namespace detail {
constexpr bool layer_catalog_well_formed = [] {
  // The 116 settings (index 0 is the sentinel), as a range projected on number.
  auto const settings = kLayerCatalog.by_neg_hall | std::views::drop(1);
  return kLayerGroupTypes.size() ==
             static_cast<std::size_t>(kNumLayerHallNumbers) + 1 &&
         // Every setting names a layer group 1..80...
         std::ranges::all_of(settings,
                             [](SpacegroupType const &t) {
                               return t.number >= 1 &&
                                      t.number <= kNumLayerGroups;
                             }) &&
         // ...grouped in non-decreasing blocks (so layer_default_hall finds the
         // first / canonical setting of each)...
         std::ranges::is_sorted(settings, {}, &SpacegroupType::number) &&
         // ...peaking at 80, with every number resolving to a default setting.
         std::ranges::max(settings, {}, &SpacegroupType::number).number ==
             kNumLayerGroups &&
         std::ranges::all_of(std::views::iota(1, kNumLayerGroups + 1),
                             [](int number) {
                               return layer_default_hall(number) != 0;
                             });
}();
} // namespace detail

static_assert(detail::layer_catalog_well_formed);
static_assert(num_layer_groups() == 80);
static_assert(layer_default_hall(1) == -1);  // p1 is the first setting
static_assert(layer_default_hall(80) != 0);  // last number resolves
static_assert(layer_default_hall(81) == 0);  // out of range
static_assert(spacegroup_type(-1).number == 1);
static_assert(spacegroup_type(-kNumLayerHallNumbers).number == kNumLayerGroups);
static_assert(spacegroup_type(0).number == 0);    // sentinel
static_assert(spacegroup_type(-117).number == 0); // out of range

[[nodiscard]] SymmetryOperations const &
operations_from_database(int hall_number);

} // namespace cppcrystal::data
