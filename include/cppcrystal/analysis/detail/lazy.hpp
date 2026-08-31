#pragma once

#include <cppcrystal/core/error.hpp>

#include <concepts>
#include <optional>
#include <type_traits>
#include <utility>

namespace cppcrystal::analysis::detail {

// A memoized Result-producing computation: `get(compute)` runs `compute` on
// the first call, caches the success value, and hands back a pointer into the
// cache (so projections copy out only the field they need). On error nothing
// is cached and the next call re-runs — boost::leaf::result is move-only, so
// the cache stores the success value, never the Result.
//
// Not race-free: concurrent first-calls race on the internal optional. The
// analyzers' warm() contract (fill every cache once, then share the const
// instance read-only) lives here in one place.
template <class T> class Lazy {
public:
  template <class F>
    requires std::same_as<std::invoke_result_t<F &>, Result<T>>
  [[nodiscard]] Result<T const *> get(F &&compute) const {
    if (!value_) {
      BOOST_LEAF_AUTO(v, compute());
      value_ = std::move(v);
    }
    return &*value_;
  }

private:
  mutable std::optional<T> value_;
};

} // namespace cppcrystal::analysis::detail
