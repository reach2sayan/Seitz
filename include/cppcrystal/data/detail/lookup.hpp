#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <functional>
#include <numeric>
#include <span>

// Shared lookup idioms of the data catalogs: every table here is 1-based with
// index 0 holding the zero-valued out-of-range fallback, and the layer-group
// family is keyed by NEGATIVE Hall numbers into a parallel table. Both
// conventions live only in the two functions below. BucketIndex is the third
// idiom: a compile-time inverted index from a small integer key to the ids
// that carry it.
namespace cppcrystal::data::detail {

// key -> the values of the ids with that key, as one contiguous span. Built
// by bucket_index at compile time; K is the number of keys (0..K-1).
template <std::size_t N, std::size_t K> struct BucketIndex {
  std::array<int, N> values{};      // grouped by key, ascending id within
  std::array<int, K + 1> offsets{}; // bucket k = values[offsets[k], offsets[k+1])

  [[nodiscard]] constexpr std::span<int const>
  operator[](int key) const noexcept {
    if (key < 0 || key >= static_cast<int>(K)) {
      return {};
    }
    auto const k = static_cast<std::size_t>(key);
    return {values.data() + offsets[k], values.data() + offsets[k + 1]};
  }
};

// Counting sort of the ids 1..N by key_of(id) into [0, K), storing
// value_of(id). Stable, so the values of one key keep ascending id order. A
// key outside [0, K) is a compile-time error.
template <std::size_t N, std::size_t K, class KeyOf,
          class ValueOf = std::identity>
[[nodiscard]] consteval BucketIndex<N, K> bucket_index(KeyOf key_of,
                                                       ValueOf value_of = {}) {
  BucketIndex<N, K> out{};
  for (int id = 1; id <= static_cast<int>(N); ++id) {
    ++out.offsets[static_cast<std::size_t>(key_of(id)) + 1];
  }
  std::partial_sum(out.offsets.begin(), out.offsets.end(),
                   out.offsets.begin());
  std::array<int, K> cursor{};
  std::copy_n(out.offsets.begin(), K, cursor.begin());
  for (int id = 1; id <= static_cast<int>(N); ++id) {
    auto const k = static_cast<std::size_t>(key_of(id));
    out.values[static_cast<std::size_t>(cursor[k]++)] =
        static_cast<int>(value_of(id));
  }
  return out;
}

// A catalog table: sized and indexable, with entry 0 as the sentinel.
template <typename T>
concept SentinelTable = requires(T const &t) {
  { t.size() } -> std::convertible_to<std::size_t>;
  t[std::size_t{0}];
};

// 1-based table lookup; index 0 is the out-of-range fallback entry.
[[nodiscard]] constexpr auto const &at_or_sentinel(SentinelTable auto const &table,
                                                   int key) noexcept {
  bool const in_range = key >= 1 && key < static_cast<int>(table.size());
  return table[static_cast<std::size_t>(in_range ? key : 0)];
}

// Dispatch a Hall number to the 3D table (positive, 1..N3) or the layer table
// (negative, keyed by -hall). The two tables must share their element type.
[[nodiscard]] constexpr auto const &hall_indexed(SentinelTable auto const &main,
                                                 SentinelTable auto const &layer,
                                                 int hall_number) noexcept {
  return hall_number < 0 ? at_or_sentinel(layer, -hall_number)
                         : at_or_sentinel(main, hall_number);
}

} // namespace cppcrystal::data::detail
