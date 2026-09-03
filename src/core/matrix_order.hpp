#pragma once

#include <cppcrystal/core/types.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#include <algorithm>
#include <concepts>
#include <functional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

// Ordering, de-duplication and lookup vocabulary for integer rotation
// matrices. Every "distinct rotations" and "find this rotation" container in
// the library is one of the sorted flat containers below, keyed by this one
// comparator: |G| <= 192, so a handful of 9-int comparisons per lookup beats
// hashing, and the ordering makes "first match" deterministic.
namespace cppcrystal {

// Strict weak order over the 9 entries of a Matrix3i (storage-order
// lexicographic; any fixed order works for set keys).
struct Matrix3iLess {
  [[nodiscard]] bool operator()(Matrix3i const &a,
                                Matrix3i const &b) const noexcept {
    return std::ranges::lexicographical_compare(std::span(a.data(), 9),
                                                std::span(b.data(), 9));
  }
};

using RotationSet = boost::container::flat_set<Matrix3i, Matrix3iLess>;
template <class T>
using RotationMap = boost::container::flat_map<Matrix3i, T, Matrix3iLess>;
template <class T>
using RotationMultimap =
    boost::container::flat_multimap<Matrix3i, T, Matrix3iLess>;

// The distinct rotations of a range. `proj` maps an element to its Matrix3i
// key (identity for ranges of rotations, &SymmetryOperation::rotation for
// operation lists). O(n log n).
template <std::ranges::input_range R, class Proj = std::identity>
[[nodiscard]] RotationSet rotation_set(R &&range, Proj proj = {}) {
  RotationSet out;
  for (auto const &item : range) {
    out.insert(std::invoke(proj, item));
  }
  return out;
}

// Whether any rotation occurs twice; stops at the first repeat.
template <std::ranges::input_range R, class Proj = std::identity>
[[nodiscard]] bool has_duplicate_rotation(R &&range, Proj proj = {}) {
  RotationSet seen;
  return std::ranges::any_of(range, [&](auto const &item) {
    return !seen.insert(std::invoke(proj, item)).second;
  });
}

// Keep the first element per distinct rotation, preserving encounter order.
template <std::ranges::input_range R, class Proj = std::identity>
[[nodiscard]] auto unique_by_rotation(R &&range, Proj proj = {})
    -> std::vector<std::ranges::range_value_t<R>> {
  std::vector<std::ranges::range_value_t<R>> out;
  RotationSet seen;
  for (auto const &item : range) {
    if (seen.insert(std::invoke(proj, item)).second) {
      out.push_back(item);
    }
  }
  return out;
}

namespace detail {
// Sort (key, index) pairs by key, then index, and adopt them as a flat
// multimap: equal_range(key) then walks the elements with that key in their
// original order, and find(key) is the first of them.
template <class Key, class Less>
[[nodiscard]] auto ordered_multimap(std::vector<std::pair<Key, int>> pairs,
                                    Less less) {
  std::ranges::sort(pairs, [&](auto const &a, auto const &b) {
    if (less(a.first, b.first)) {
      return true;
    }
    return !less(b.first, a.first) && a.second < b.second;
  });
  return boost::container::flat_multimap<Key, int, Less>(
      boost::container::ordered_range, pairs.begin(), pairs.end());
}
} // namespace detail

// rotation -> indices of the elements carrying it, ascending within a
// rotation. O(n log n) to build.
template <std::ranges::random_access_range R, class Proj = std::identity>
[[nodiscard]] RotationMultimap<int> index_by_rotation(R const &range,
                                                      Proj proj = {}) {
  std::vector<std::pair<Matrix3i, int>> pairs;
  pairs.reserve(std::ranges::size(range));
  for (auto const [i, item] : range | std::views::enumerate) {
    pairs.emplace_back(std::invoke(proj, item), static_cast<int>(i));
  }
  return detail::ordered_multimap(std::move(pairs), Matrix3iLess{});
}

// The exactly-compared part of a magnetic operation: rotation plus
// time-reversal flag. Translations are compared with tolerance inside the
// bucket an OperationKey selects.
struct OperationKey {
  Matrix3i rotation;
  bool time_reversal;
};

struct OperationKeyLess {
  [[nodiscard]] bool operator()(OperationKey const &a,
                                OperationKey const &b) const noexcept {
    if (a.time_reversal != b.time_reversal) {
      return !a.time_reversal;
    }
    return Matrix3iLess{}(a.rotation, b.rotation);
  }
};

template <class Op>
concept MagneticOperationLike = requires(Op const &op) {
  { op.spatial.rotation } -> std::convertible_to<Matrix3i>;
  { op.time_reversal } -> std::convertible_to<bool>;
};

using OperationMultimap =
    boost::container::flat_multimap<OperationKey, int, OperationKeyLess>;

// (rotation, time_reversal) -> indices, ascending within a key.
template <std::ranges::random_access_range R>
  requires MagneticOperationLike<std::ranges::range_value_t<R>>
[[nodiscard]] OperationMultimap index_by_operation_key(R const &range) {
  std::vector<std::pair<OperationKey, int>> pairs;
  pairs.reserve(std::ranges::size(range));
  for (auto const [i, op] : range | std::views::enumerate) {
    pairs.emplace_back(OperationKey{op.spatial.rotation, op.time_reversal},
                       static_cast<int>(i));
  }
  return detail::ordered_multimap(std::move(pairs), OperationKeyLess{});
}

// Append `value` unless an element equivalent under `equiv` is already
// present; returns whether it was appended. The tolerance-aware analogue of
// ranges::unique for an unsorted container — first occurrence wins.
template <std::ranges::forward_range C, class Equiv>
bool push_unique(C &out, std::ranges::range_value_t<C> value, Equiv &&equiv) {
  if (std::ranges::any_of(out,
                          [&](auto const &e) { return equiv(e, value); })) {
    return false;
  }
  out.push_back(std::move(value));
  return true;
}

} // namespace cppcrystal
