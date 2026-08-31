#pragma once

#include <cppcrystal/core/types.hpp>

#include <boost/container/flat_set.hpp>

#include <algorithm>
#include <functional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

// Ordering and de-duplication vocabulary for integer rotation matrices. Every
// "distinct rotations" container in the library shares this one comparator.
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

// Keep the first element per distinct rotation, preserving encounter order.
// `proj` maps an element to its Matrix3i key (identity for ranges of
// rotations, &SymmetryOperation::rotation for operation lists). O(n log n).
template <std::ranges::input_range R, class Proj = std::identity>
[[nodiscard]] auto unique_by_rotation(R &&range, Proj proj = {})
    -> std::vector<std::ranges::range_value_t<R>> {
  std::vector<std::ranges::range_value_t<R>> out;
  boost::container::flat_set<Matrix3i, Matrix3iLess> seen;
  for (auto const &item : range) {
    if (seen.insert(std::invoke(proj, item)).second) {
      out.push_back(item);
    }
  }
  return out;
}

// Append `value` unless an element equivalent under `equiv` is already
// present; returns whether it was appended. The tolerance-aware analogue of
// ranges::unique for an unsorted container — first occurrence wins.
template <class T, class Equiv>
bool push_unique(std::vector<T> &out, T value, Equiv &&equiv) {
  if (std::ranges::any_of(out, [&](T const &e) { return equiv(e, value); })) {
    return false;
  }
  out.push_back(std::move(value));
  return true;
}

} // namespace cppcrystal
