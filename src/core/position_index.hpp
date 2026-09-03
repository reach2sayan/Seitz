#pragma once

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

#include <boost/geometry/algorithms/disjoint.hpp>
#include <boost/geometry/algorithms/intersects.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <concepts>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

namespace cppcrystal {

// True when fractional points a and b are the same site: their difference,
// folded to the minimal image along the periodic axes only, is within a
// Cartesian distance of symprec (inclusive). The one definition of "same site"
// for any periodicity.
[[nodiscard]] bool coincident(Vector3d const &a, Vector3d const &b,
                              Matrix3d const &lattice, double symprec,
                              CellPeriodicity const &periodicity) noexcept;

// A build-once, query-many index answering "which atoms sit at this
// fractional point". An R-tree over the atoms' Cartesian positions, each
// folded into the cell along the periodic axes. A query folds its point the
// same way and asks the tree for everything inside a symprec-sized box around
// each image the minimal-image fold in `coincident` can pick: shifts of -1, 0,
// +1 along every periodic axis. The symprec sphere sits inside that box, so
// candidates() is a guaranteed superset of the coincident atoms; matches()
// re-tests each with the exact predicate, so a caller replaces its linear
// scan without changing what it accepts. Candidates come out ascending, so a
// first hit is the first hit of the linear scan it replaces.
//
// Non-owning: the positions and types must outlive the index.
class PositionIndex {
public:
  PositionIndex(Positions const &positions, Types const &types,
                Matrix3d const &lattice, double symprec,
                CellPeriodicity const &periodicity);
  PositionIndex(Cell const &cell, double symprec);

  PositionIndex(PositionIndex const &) = delete;
  PositionIndex &operator=(PositionIndex const &) = delete;

  // Atoms in the box neighbourhood of `point`: a superset of the coincident
  // ones, ascending, each once.
  [[nodiscard]] std::vector<int> candidates(Vector3d const &point) const;

  // Atoms coincident with `point`, ascending.
  [[nodiscard]] auto matches(Vector3d const &point) const {
    return candidates(point) | std::views::filter([this, point](int atom) {
             return coincides(point, atom);
           });
  }

  // Atoms of `type` coincident with `point`, ascending.
  [[nodiscard]] auto matches(Vector3d const &point, int type) const {
    return candidates(point) |
           std::views::filter([this, point, type](int atom) {
             return type_of(atom) == type && coincides(point, atom);
           });
  }

  // The lowest-index atom of `type` coincident with `point` that `accept`s:
  // the first-hit-in-index-order primitive of the linear scans this replaces.
  [[nodiscard]] std::optional<int>
  first_match(Vector3d const &point, int type,
              std::predicate<int> auto accept) const {
    for (int const atom : matches(point, type)) {
      if (accept(atom)) {
        return atom;
      }
    }
    return std::nullopt;
  }
  [[nodiscard]] std::optional<int> first_match(Vector3d const &point,
                                               int type) const {
    return first_match(point, type, [](int) { return true; });
  }

private:
  using Point =
      boost::geometry::model::point<double, 3, boost::geometry::cs::cartesian>;
  using Box = boost::geometry::model::box<Point>;
  using Tree = boost::geometry::index::rtree<
      std::pair<Point, int>, boost::geometry::index::quadratic<16>>;

  // Cartesian position of a fractional point folded into the cell.
  [[nodiscard]] Vector3d cartesian(Vector3d const &frac) const noexcept;
  [[nodiscard]] bool coincides(Vector3d const &point, int atom) const noexcept;
  [[nodiscard]] int type_of(int atom) const noexcept {
    return (*types_)[static_cast<std::size_t>(atom)];
  }

  Matrix3d lattice_;
  double symprec_;
  double half_width_; // of the query box: symprec plus a rounding margin
  CellPeriodicity periodicity_;
  Positions const *positions_;
  Types const *types_;
  Tree tree_;
};

} // namespace cppcrystal
