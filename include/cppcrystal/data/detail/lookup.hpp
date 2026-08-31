#pragma once

#include <cstddef>

// Shared lookup idioms of the data catalogs: every table here is 1-based with
// index 0 holding the zero-valued out-of-range fallback, and the layer-group
// family is keyed by NEGATIVE Hall numbers into a parallel table. Both
// conventions live only in these two functions.
namespace cppcrystal::data::detail {

// 1-based table lookup; index 0 is the out-of-range fallback entry.
template <typename Table>
[[nodiscard]] constexpr auto const &at_or_sentinel(Table const &table,
                                                   int key) noexcept {
  bool const in_range = key >= 1 && key < static_cast<int>(table.size());
  return table[static_cast<std::size_t>(in_range ? key : 0)];
}

// Dispatch a Hall number to the 3D table (positive, 1..N3) or the layer table
// (negative, keyed by -hall). The two tables must share their element type.
template <typename Main, typename Layer>
[[nodiscard]] constexpr auto const &
hall_indexed(Main const &main, Layer const &layer, int hall_number) noexcept {
  return hall_number < 0 ? at_or_sentinel(layer, -hall_number)
                         : at_or_sentinel(main, hall_number);
}

} // namespace cppcrystal::data::detail
