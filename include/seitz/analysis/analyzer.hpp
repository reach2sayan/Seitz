#pragma once

#include <seitz/core/error.hpp>

#include <atomic>
#include <concepts>
#include <mutex>
#include <optional>
#include <type_traits>
#include <utility>

#pragma GCC visibility push(default)

namespace seitz::analysis {

namespace detail {

// A memoized Result-producing computation:
// Race-free: a populated cache is read lock-free through the acquire flag;
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
template <typename T>
concept AnalyzerTraits = requires {
  typename T::CellType;
  typename T::DatasetType;
  typename T::ToleranceType;
};

// A persistent, stateful view over a cell + tolerances that lazily computes
// and memoizes its determination.
// Facade over the pipeline, Template Method over the determination:
// owns the inputs, the cache and the projection
// machinery; `Derived` supplies only `determine()` (the actual pipeline).
// Thread-safety: every memo is race-free (detail::Lazy)
template <typename Derived, AnalyzerTraits Traits> class Analyzer {
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
  // rather than a dangling read.
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

} // namespace seitz::analysis

#pragma GCC visibility pop
