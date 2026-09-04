#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <numeric>
#include <span>

// The inverted index the data catalogs are queried through: a compile-time map
// from a small integer key to the values that carry it, built by counting sort
// at consteval time. (The 1-based-with-sentinel-row and negative-Hall lookups
// that used to live here are gone: Catalog<Family> is keyed by a validated
// HallNumber, so every lookup is in range by construction.)
namespace cppcrystal::data::detail {

// key -> the values of the ids with that key, as one contiguous span. Built by
// bucket_index at compile time; K is the number of keys (0..K-1). T is usually
// an id, but any literal type works -- the t-subgroup adjacency in
// group/subgroup_graph.hpp buckets whole relation records this way.
template <class T, std::size_t N, std::size_t K> struct BucketIndex {
  std::array<T, N> values{}; // grouped by key, ascending id within
  std::array<int, K + 1>
      offsets{}; // bucket k = values[offsets[k], offsets[k+1])

  [[nodiscard]] constexpr std::span<T const>
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
template <std::size_t N, std::size_t K, class T = int, class KeyOf,
          class ValueOf = std::identity>
[[nodiscard]] consteval BucketIndex<T, N, K>
bucket_index(KeyOf key_of, ValueOf value_of = {}) {
  BucketIndex<T, N, K> out{};
  for (int id = 1; id <= static_cast<int>(N); ++id) {
    ++out.offsets[static_cast<std::size_t>(key_of(id)) + 1];
  }
  std::partial_sum(out.offsets.begin(), out.offsets.end(), out.offsets.begin());
  std::array<int, K> cursor{};
  std::copy_n(out.offsets.begin(), K, cursor.begin());
  for (int id = 1; id <= static_cast<int>(N); ++id) {
    auto const k = static_cast<std::size_t>(key_of(id));
    out.values[static_cast<std::size_t>(cursor[k]++)] = value_of(id);
  }
  return out;
}

} // namespace cppcrystal::data::detail
