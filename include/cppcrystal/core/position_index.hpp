#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

#include <boost/container/static_vector.hpp>

#include <algorithm>
#include <array>
#include <compare>
#include <concepts>
#include <cstdint>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

namespace cppcrystal {

// True when fractional points a and b are the same site: their difference,
// folded to the minimal image along the periodic axes only, is within a
// Cartesian distance of symprec (inclusive). The one definition of "same site"
// for any periodicity; is_overlap is its single-aperiodic-axis spelling.
[[nodiscard]] bool coincident(Vector3d const &a, Vector3d const &b,
                              Matrix3d const &lattice, double symprec,
                              CellPeriodicity const &periodicity) noexcept;

// Bucket grid over fractional space, sized so that two points within symprec
// (Cartesian) land in the same or an adjacent bucket along every axis:
//   |Δfrac_i| = |row_i(L⁻¹) · Δcart| ≤ ‖row_i(L⁻¹)‖ · symprec,
// and a bucket edge is at least twice that bound, so float rounding at a
// bucket boundary can never push a match two buckets away. Periodic axes are
// cut into `divisions` buckets that wrap; an aperiodic axis is an unbounded
// grid of `width`-sized buckets with no wrap.
struct BucketGeometry {
  using Key = std::array<std::int64_t, 3>;

  Vector3d width;      // fractional edge per axis
  Vector3i divisions;  // >= 1 on a periodic axis; 0 marks an aperiodic axis

  [[nodiscard]] static BucketGeometry
  of(Matrix3d const &lattice, double symprec,
     CellPeriodicity const &periodicity) noexcept;

  [[nodiscard]] Key bucket_of(Vector3d const &frac) const noexcept;
};

// A build-once, query-many index answering "which atoms sit at this
// fractional point". Sorted (bucket, atom) pairs queried by equal_range over
// the 3^3 bucket neighbourhood: candidates() is a guaranteed superset of the
// coincident atoms, matches() re-tests each with the exact `coincident`
// predicate, so a caller replaces its linear scan without changing what it
// accepts. Ascending atom index within a bucket keeps first-match queries
// deterministic. O(n log n) to build, O(27 log n + hits) per query.
//
// Non-owning: the positions and types must outlive the index.
class PositionIndex {
public:
  PositionIndex(BucketGeometry geometry, Positions const &positions,
                Types const &types, Matrix3d const &lattice, double symprec,
                CellPeriodicity const &periodicity);
  PositionIndex(Cell const &cell, double symprec);

  PositionIndex(PositionIndex const &) = delete;
  PositionIndex &operator=(PositionIndex const &) = delete;

  // Atoms in the bucket neighbourhood of `point`: a superset of the
  // coincident ones, in no particular order.
  [[nodiscard]] auto candidates(Vector3d const &point) const {
    return bucket_ranges(point) | std::views::join |
           std::views::transform(&Entry::atom);
  }

  // Atoms coincident with `point`, in no particular order.
  [[nodiscard]] auto matches(Vector3d const &point) const {
    return candidates(point) | std::views::filter([this, point](int atom) {
             return coincides(point, atom);
           });
  }

  // Atoms of `type` coincident with `point`, in no particular order.
  [[nodiscard]] auto matches(Vector3d const &point, int type) const {
    return candidates(point) | std::views::filter([this, point, type](int atom) {
             return type_of(atom) == type && coincides(point, atom);
           });
  }

  // The lowest-index atom of `type` coincident with `point` that `accept`s:
  // the first-hit-in-index-order primitive of the linear scans this replaces.
  [[nodiscard]] std::optional<int>
  first_match(Vector3d const &point, int type,
              std::predicate<int> auto accept) const {
    std::optional<int> best;
    for (int const atom : matches(point, type)) {
      if (accept(atom) && (!best || atom < *best)) {
        best = atom;
      }
    }
    return best;
  }
  [[nodiscard]] std::optional<int> first_match(Vector3d const &point,
                                               int type) const {
    return first_match(point, type, [](int) { return true; });
  }

private:
  struct Entry {
    BucketGeometry::Key bucket;
    int atom;
    auto operator<=>(Entry const &) const = default;
  };
  using Buckets = boost::container::static_vector<std::span<Entry const>, 27>;

  // The non-empty entry runs of the 3^3 neighbourhood, each bucket once.
  [[nodiscard]] Buckets bucket_ranges(Vector3d const &point) const;
  [[nodiscard]] bool coincides(Vector3d const &point,
                               int atom) const noexcept;
  [[nodiscard]] int type_of(int atom) const noexcept {
    return (*types_)[static_cast<std::size_t>(atom)];
  }

  BucketGeometry geometry_;
  Matrix3d lattice_;
  double symprec_;
  CellPeriodicity periodicity_;
  Positions const *positions_;
  Types const *types_;
  std::vector<Entry> entries_; // sorted by (bucket, atom)
};

} // namespace cppcrystal
