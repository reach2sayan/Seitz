#pragma once

#include <spglib/core/symmetry_operation.hpp>
#include <spglib/data/rod_group_tables.hpp>

#include <cstddef>
#include <ranges>
#include <string_view>

namespace spglib::data {

// Whether `rod_number` is a valid rod-group number (1..75).
[[nodiscard]] constexpr bool rod_number_in_range(int rod_number) noexcept {
  return rod_number >= 1 && rod_number <= kNumRodGroups;
}

[[nodiscard]] constexpr int num_rod_groups() noexcept { return kNumRodGroups; }
[[nodiscard]] constexpr int rod_periodic_axis() noexcept {
  return kRodPeriodicAxis;
}

// Hermann-Mauguin symbol of rod group `rod_number` (1..75); empty out of range.
[[nodiscard]] constexpr std::string_view rod_symbol(int rod_number) noexcept {
  if (!rod_number_in_range(rod_number)) {
    return {};
  }
  return kRodGroupSymbols[static_cast<std::size_t>(rod_number - 1)];
}

namespace detail {
constexpr bool rod_offsets_well_formed = [] {
  if (kRodGroupSymbols.size() != static_cast<std::size_t>(kNumRodGroups)) {
    return false;
  }
  if (kRodOperationOffset.front() != 0 ||
      kRodOperationOffset.back() != static_cast<int>(kRodOperations.size())) {
    return false;
  }
  return std::ranges::adjacent_find(kRodOperationOffset,
                                    std::ranges::greater_equal{}) ==
         kRodOperationOffset.end();
}();
} // namespace detail

static_assert(detail::rod_offsets_well_formed);
static_assert(num_rod_groups() == 75);
static_assert(rod_periodic_axis() == 2);
static_assert(rod_symbol(0).empty());  // out of range (below)
static_assert(rod_symbol(76).empty()); // out of range (above)
static_assert(rod_symbol(1) == "p1");
static_assert(rod_symbol(75) == "p6/mmc");

[[nodiscard]] SymmetryOperations rod_operations_from_database(int rod_number);

} // namespace spglib::data
