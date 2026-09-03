#pragma once

#include <cppcrystal/core/error.hpp>

#include <concepts>
#include <optional>
#include <type_traits>
#include <utility>

namespace cppcrystal::analysis {

namespace detail {

// A memoized Result-producing computation: `get(compute)` runs `compute` on
// the first call, caches the success value, and hands back a pointer into the
// cache (so projections copy out only the field they need). On error nothing
// is cached and the next call re-runs — boost::leaf::result is move-only, so
// the cache stores the success value, never the Result.
//
// Not race-free: concurrent first-calls race on the internal optional. The
// analyzers' warm() contract (fill every cache once, then share the const
// instance read-only) is what makes that safe.
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

} // namespace detail

// A persistent, stateful view over a cell + tolerances that lazily computes
// and memoizes its determination. Facade over the pipeline, Template Method
// over the determination: this owns the inputs, the cache and the projection
// machinery; `Derived` supplies only `determine()` (the actual pipeline) and
// `warm_derived()` (its own extra caches, if any).
//
// Static polymorphism only — Derived is a concrete type, there is no runtime
// hierarchy and no virtual dispatch.
//
// Thread-safety: the per-instance caches are NOT race-free. To share one
// instance read-only across threads, call warm() once on a single thread
// first; afterwards every getter is served from a populated cache and does no
// writing. (The shared global tables are primed by cppcrystal::warmup().)
template <class Derived, class Traits> class Analyzer {
public:
  using CellType = typename Traits::CellType;
  using DatasetType = typename Traits::DatasetType;
  using ToleranceType = typename Traits::ToleranceType;

  [[nodiscard]] CellType const &cell() const noexcept { return cell_; }
  [[nodiscard]] ToleranceType const &tolerance() const noexcept {
    return tol_;
  }

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

  // Force every lazy cache. Returns the first error encountered, or success
  // once all caches are populated; afterwards this const instance may be
  // shared read-only across threads.
  Result<void> warm() const {
    BOOST_LEAF_CHECK(cached_dataset());
    return derived().warm_derived();
  }

protected:
  Analyzer(CellType cell, ToleranceType tol)
      : cell_(std::move(cell)), tol_(std::move(tol)) {}

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

  // Nothing extra to warm unless Derived says otherwise.
  [[nodiscard]] Result<void> warm_derived() const { return {}; }

  CellType cell_;
  ToleranceType tol_;

private:
  detail::Lazy<DatasetType> dataset_;
};

} // namespace cppcrystal::analysis
