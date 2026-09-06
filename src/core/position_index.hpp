#pragma once

#include "core/testable.hpp"
#include <seitz/core/cell.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>

#include <boost/container/small_vector.hpp>

#include <array>
#include <concepts>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

namespace seitz {

// True when fractional points a and b are the same site: their difference,
// folded to the minimal image along the periodic axes only, is within a
// Cartesian distance of symprec (inclusive). The one definition of "same site"
// for any periodicity.
// Inline: it is the innermost test of every search, and out of line its 3x3
// product and fold were half again the cost of the call itself.
[[nodiscard]] inline bool coincident(Vector3d const &a, Vector3d const &b,
                                     Matrix3d const &lattice, double symprec,
                                     CellPeriodicity const &periodicity) noexcept {
  return (lattice * minimal_image(a - b, periodicity)).norm() <= symprec;
}

// Build-once, query-many index of "which atoms sit at this fractional point":
// a uniform grid over the cell, atoms bucketed by the cell their folded
// fractional position falls in, about one atom per bucket. A query folds the
// same way and scans every bucket that `coincident`'s minimal-image fold could
// land a match in: a Cartesian offset of at most half_width is a fractional
// one of at most ||lattice^-1||_inf * half_width (image_reach_), so those are
// the buckets within image_reach_ of the point along each periodic axis,
// wrapping at the cell boundary -- typically one, never more than eight for a
// sane lattice. The buckets over-approximate the sphere, so candidates() is a
// superset and matches() re-tests exactly; buckets hold atoms in index order
// and a multi-bucket scan is sorted, so candidates come out ascending and a
// first hit is the linear scan's first hit. An aperiodic axis has a single
// bucket: its coordinate is unbounded and unfolded.
//
// This replaced a Boost.Geometry R-tree once the algorithmic work made the
// query itself the hot spot: the tree cost ~150 ns per lookup and a sizeable
// build even for 16 atoms, the grid a few tens of ns and a counting sort.
//
// Queries take a caller-owned Scratch: the hot paths run once per atom per
// candidate operation and must not allocate.
//
// Non-owning: the positions and types must outlive the index.
class SEITZ_TESTABLE PositionIndex {
public:
  // Reused across queries; each query clears it. Sized for the common case:
  // a coincidence query typically returns one or two atoms.
  using Scratch = boost::container::small_vector<int, 8>;

  PositionIndex(Positions const &positions, Types const &types,
                Matrix3d const &lattice, double symprec,
                CellPeriodicity const &periodicity);
  PositionIndex(Cell const &cell, double symprec);

  PositionIndex(PositionIndex const &) = delete;
  PositionIndex &operator=(PositionIndex const &) = delete;

  // Atoms in the bucket neighbourhood of `point`: a superset of the coincident
  // ones, ascending, each once. A view of `out`, which it refills -- or of
  // the index's own all-atoms list when the neighbourhood is the whole cell.
  [[nodiscard]] std::span<int const> candidates(Vector3d const &point,
                                                Scratch &out) const;

  // Allocating form, for callers not in a loop.
  [[nodiscard]] std::vector<int> candidates(Vector3d const &point) const;

  // Atoms coincident with `point`, ascending. The view borrows `out`, which
  // must outlive it and must not be reused while it is being iterated.
  [[nodiscard]] auto matches(Vector3d const &point, Scratch &out) const {
    return candidates(point, out) | std::views::filter([this, point](int atom) {
             return coincides(point, atom);
           });
  }
  [[nodiscard]] auto matches(Vector3d const &point) const {
    return candidates(point) | std::views::filter([this, point](int atom) {
             return coincides(point, atom);
           });
  }

  // Atoms of `type` coincident with `point`, ascending.
  [[nodiscard]] auto matches(Vector3d const &point, int type,
                             Scratch &out) const {
    return candidates(point, out) |
           std::views::filter([this, point, type](int atom) {
             return type_of(atom) == type && coincides(point, atom);
           });
  }
  [[nodiscard]] auto matches(Vector3d const &point, int type) const {
    return candidates(point) |
           std::views::filter([this, point, type](int atom) {
             return type_of(atom) == type && coincides(point, atom);
           });
  }

  // The lowest-index atom of `type` coincident with `point` that `accept`s:
  // the first-hit-in-index-order primitive of the linear scans this replaces.
  [[nodiscard]] std::optional<int>
  first_match(Vector3d const &point, int type, Scratch &out,
              std::predicate<int> auto accept) const {
    for (int const atom : matches(point, type, out)) {
      if (accept(atom)) {
        return atom;
      }
    }
    return std::nullopt;
  }
  [[nodiscard]] std::optional<int> first_match(Vector3d const &point, int type,
                                               Scratch &out) const {
    return first_match(point, type, out, [](int) { return true; });
  }
  [[nodiscard]] std::optional<int>
  first_match(Vector3d const &point, int type,
              std::predicate<int> auto accept) const {
    Scratch out;
    return first_match(point, type, out, accept);
  }
  [[nodiscard]] std::optional<int> first_match(Vector3d const &point,
                                               int type) const {
    return first_match(point, type, [](int) { return true; });
  }

private:
  // Bucket along `axis` of a folded fractional coordinate, wrapped into
  // [0, divisions_[axis]).
  [[nodiscard]] int bucket_along(double x, std::size_t axis) const noexcept;
  [[nodiscard]] bool coincides(Vector3d const &point, int atom) const noexcept {
    return coincident(point, positions_->row(atom).transpose(), lattice_,
                      symprec_, periodicity_);
  }
  [[nodiscard]] int type_of(int atom) const noexcept {
    return (*types_)[static_cast<std::size_t>(atom)];
  }

  Matrix3d lattice_;
  double symprec_;
  double image_reach_; // fractional reach of a coincident atom, per axis
  CellPeriodicity periodicity_;
  Positions const *positions_;
  Types const *types_;
  std::array<int, 3> divisions_; // buckets per axis; 1 on an aperiodic axis
  std::vector<int> starts_;      // bucket -> first entry of atoms_, plus end
  std::vector<int> atoms_;       // atom indices by bucket, ascending within
  std::vector<int> everything_;  // 0..n-1, the answer when every bucket is hit
};

} // namespace seitz
