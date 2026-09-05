#pragma once

#include <cppcrystal/core/error.hpp>

#include <atomic>
#include <concepts>
#include <mutex>
#include <optional>
#include <type_traits>
#include <utility>

#pragma GCC visibility push(default)

namespace cppcrystal::analysis {

namespace detail {

// A memoized Result-producing computation: `get(compute)` runs `compute` on
// the first call, caches the success value, and hands back a pointer into the
// cache. On error nothing is cached and the next call re-runs —
// boost::leaf::result is move-only, so the cache stores the success value,
// never the Result.
//
// Race-free: a populated cache is read lock-free through the acquire flag;
// first calls serialise on the mutex, so concurrent callers of a shared const
// analyzer see exactly one computation. Moving a Lazy moves the value and
// starts a fresh guard (moving while another thread reads is a caller bug,
// as for any object).
template <class T> class Lazy {
public:
  Lazy() = default;
  Lazy(Lazy &&other) noexcept
      : value_(std::move(other.value_)), ready_(value_.has_value()) {}
  Lazy &operator=(Lazy &&other) noexcept {
    value_ = std::move(other.value_);
    ready_.store(value_.has_value(), std::memory_order_release);
    return *this;
  }
  Lazy(Lazy const &) = delete;
  Lazy &operator=(Lazy const &) = delete;

  template <class F>
    requires std::same_as<std::invoke_result_t<F &>, Result<T>>
  [[nodiscard]] Result<T const *> get(F &&compute) const {
    if (ready_.load(std::memory_order_acquire)) {
      return &*value_;
    }
    std::lock_guard const lock(mutex_);
    if (!value_) {
      BOOST_LEAF_AUTO(v, compute());
      value_ = std::move(v);
      ready_.store(true, std::memory_order_release);
    }
    return &*value_;
  }

private:
  mutable std::optional<T> value_;
  mutable std::atomic<bool> ready_{false};
  mutable std::mutex mutex_;
};

} // namespace detail

// What an analyzer's Traits names: the cell it takes, the dataset its
// determination produces, and the tolerance that drives it.
template <class T>
concept AnalyzerTraits = requires {
  typename T::CellType;
  typename T::DatasetType;
  typename T::ToleranceType;
};

// A persistent, stateful view over a cell + tolerances that lazily computes
// and memoizes its determination. Facade over the pipeline, Template Method
// over the determination: this owns the inputs, the cache and the projection
// machinery; `Derived` supplies only `determine()` (the actual pipeline).
//
// Static polymorphism only — Derived is a concrete type, there is no runtime
// hierarchy and no virtual dispatch.
//
// Thread-safety: every memo is race-free (detail::Lazy), so a const analyzer
// may be shared across threads from the moment it is built; the first caller
// of each query pays for it, the rest read the cache.
template <class Derived, AnalyzerTraits Traits> class Analyzer {
public:
  using CellType = typename Traits::CellType;
  using DatasetType = typename Traits::DatasetType;
  using ToleranceType = typename Traits::ToleranceType;

  [[nodiscard]] CellType const &cell() const noexcept { return cell_; }
  [[nodiscard]] ToleranceType const &tolerance() const noexcept { return tol_; }

  // The full determination of the input cell (memoized). A reference into the
  // memo: every projection below is a view onto this one computation, and a
  // caller that wants ownership copies for itself.
  //
  // Ref-qualified, so a reference into a temporary analyzer is a compile error
  // rather than a dangling read. An analyzer exists to be held anyway --
  // `from_cell(cell).dataset()` throws the memo away with the temporary.
  [[nodiscard]] Result<DatasetType const &> dataset() const & {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return *ds;
  }
  Result<DatasetType const &> dataset() const && = delete;

protected:
  Analyzer(CellType cell, ToleranceType tol)
      : cell_{std::move(cell)}, tol_{std::move(tol)} {}

  [[nodiscard]] Derived const &derived() const noexcept {
    return static_cast<Derived const &>(*this);
  }

  [[nodiscard]] Result<DatasetType const *> cached_dataset() const {
    return dataset_.get([&] { return derived().determine(); });
  }

  // One field of the memoized dataset, by reference. Same lifetime contract as
  // dataset(): the derived accessors that use this are ref-qualified too.
  template <auto Member>
  [[nodiscard]] auto project() const
      -> Result<decltype(std::declval<DatasetType const &>().*Member) const &> {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return ds->*Member;
  }

  CellType cell_;
  ToleranceType tol_;

private:
  detail::Lazy<DatasetType> dataset_;
};

} // namespace cppcrystal::analysis

#pragma GCC visibility pop
