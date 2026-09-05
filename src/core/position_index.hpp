#pragma once

#include "core/testable.hpp"
#include <seitz/core/cell.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>

#include <boost/container/small_vector.hpp>
#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

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
[[nodiscard]] SEITZ_TESTABLE bool
coincident(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
           double symprec, CellPeriodicity const &periodicity) noexcept;

// A build-once, query-many index answering "which atoms sit at this
// fractional point". An R-tree over the atoms' Cartesian positions, each
// folded into the cell along the periodic axes. A query folds its point the
// same way and asks the tree for everything inside a symprec-sized box around
// each image the minimal-image fold in `coincident` can pick. The symprec
// sphere sits inside that box, so candidates() is a guaranteed superset of the
// coincident atoms; matches() re-tests each with the exact predicate, so a
// caller replaces its linear scan without changing what it accepts. Candidates
// come out ascending, so a first hit is the first hit of the linear scan it
// replaces.
//
// Only the images that can actually contain a match are queried. Both the
// query point and the indexed atoms are folded into [0, 1), and a Cartesian
// offset of at most `half_width_` is a fractional offset of at most
// ||lattice^-1||_inf * half_width_ -- so the +1 image along an axis can only
// match when the folded coordinate is that close to 0, and the -1 image only
// when it is that close to 1. For a typical symprec that is one box per query
// rather than the 27 of the full +-1 cube.
//
// Queries take a caller-owned Scratch buffer so the hot paths -- which query
// once per atom per candidate operation -- do not allocate. The allocating
// overloads remain for the cold callers.
//
// A uniform bucket grid was tried here and REJECTED on measurement: with the
// query box far smaller than a bucket, addressing the bucket arithmetically
// beat the tree walk by ~5% on the determination driver, but cost ~14% on
// generation, where the index is rebuilt per attempt over a small cell and the
// bucket-array build never amortizes. The tree wins on the mix; the node size
// below is the part of that experiment worth keeping.
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

  // Atoms in the box neighbourhood of `point`: a superset of the coincident
  // ones, ascending, each once. Refills `out` and returns a view of it.
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
  using Point =
      boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
  using Box = boost::geometry::model::box<Point>;
  // Node capacity 8, measured rather than assumed: against the default 16 it
  // is ~1% on the determination driver and ~12% on generation, where the cells
  // are small (12-60 atoms) and a fatter node costs more box tests per level
  // than the shallower tree saves. 4 and 32 were both worse (32 by 14%), and
  // rstar / linear were indistinguishable from quadratic at 16.
  using Tree =
      boost::geometry::index::rtree<std::pair<Point, int>,
                                    boost::geometry::index::quadratic<8>>;

  // Cartesian position of a fractional point folded into the cell.
  [[nodiscard]] Vector3d cartesian(Vector3d const &frac) const noexcept;
  [[nodiscard]] bool coincides(Vector3d const &point, int atom) const noexcept;
  [[nodiscard]] int type_of(int atom) const noexcept {
    return (*types_)[static_cast<std::size_t>(atom)];
  }

  Matrix3d lattice_;
  Matrix3d lattice_inv_;
  double symprec_;
  double half_width_;  // of the query box: symprec plus a rounding margin
  double image_reach_; // half_width_ expressed as a fractional coordinate
  CellPeriodicity periodicity_;
  Positions const *positions_;
  Types const *types_;
  Tree tree_;
};

} // namespace seitz
