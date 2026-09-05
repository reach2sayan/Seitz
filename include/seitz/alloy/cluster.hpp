#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <functional>
#include <ranges>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::alloy {

// Coincidence tolerance for cluster points, in fractional coordinates.
//
// Cluster geometry is produced by applying exact symmetry operations to exact
// lattice positions, so a true coincidence is exact to within rounding and the
// only residual of any size is the 1/n averaged translation correction in
// site_symmetry (order 1e-16). This value is therefore slack rather than a
// threshold anything lands near; it matches ATAT's global tolerance so the
// published cluster multiplicities are reproduced exactly.
inline constexpr double kClusterPrec = 1e-3;

// One point of a cluster: where it sits, how many species its sublattice
// admits, and which point function decorates it.
//
// `species` is the count itself, never ATAT's biased "site type" -- see the
// note on SiteBasis. `function` indexes the site basis, in [0, species - 1);
// a species assignment is NOT stored here but in Configuration::occupation,
// because ATAT overloading one field with both meanings is what made its
// v-matrix builder read two different quantities off the same name.
struct ClusterPoint {
  Vector3d position{Vector3d::Zero()}; // fractional, in the parent's frame
  int species = 2;
  int function = 0;
};

// A decorated cluster: a set of points, each carrying its point function.
//
// One container of points, not the three parallel vectors ATAT keeps in step
// by hand -- which is why this header needs no proxy reference, no proxy
// iterator, no tuple protocol and no iterator-interface library to make
// `for (auto const &p : cluster)` work.
//
// small_vector<_, 4>: the nearest-neighbour tetrahedron is the workhorse
// maximal cluster of the CVM, and the enumeration builds millions of these, so
// four points inline keeps the allocator out of the inner loops.
class Cluster {
public:
  using Points = boost::container::small_vector<ClusterPoint, 4>;

  Cluster() = default;
  explicit Cluster(Points points) noexcept : points_{std::move(points)} {}

  [[nodiscard]] auto begin() const noexcept { return points_.begin(); }
  [[nodiscard]] auto end() const noexcept { return points_.end(); }
  [[nodiscard]] int size() const noexcept {
    return static_cast<int>(points_.size());
  }
  [[nodiscard]] bool empty() const noexcept { return points_.empty(); }
  [[nodiscard]] ClusterPoint const &operator[](int i) const noexcept {
    return points_[static_cast<std::size_t>(i)];
  }

  // Every point moved by `op`; species and function ride along in order.
  [[nodiscard]] Cluster transformed(SymmetryOperation const &op) const;

  // Longest Cartesian distance between two of the points -- the cluster
  // "length" that the radius cutoffs and the subcluster ordering compare.
  [[nodiscard]] double diameter(Lattice const &lattice) const;

  // The same decorated cluster modulo a lattice translation: some whole-cluster
  // shift makes every point of *this coincide with a point of `other` carrying
  // the same function. In the parent's fractional frame "modulo the cell" is
  // "modulo 1", so this is a rounding test rather than a metric one.
  //
  // `species` is deliberately NOT compared, matching ATAT: a decorated cluster
  // may legitimately match one whose points sit on differently-sized
  // sublattices, and the function index is what identifies the basis element.
  [[nodiscard]] bool congruent(Cluster const &other, double tol) const;

  // The cluster as a Cell, so it can be viewed, exported, or handed to
  // SymmetryAnalyzer like any other structure in this library. Atom types are
  // the point FUNCTIONS, which makes the determination of the resulting cell
  // the symmetry of the decorated cluster rather than of its bare geometry.
  [[nodiscard]] Cell as_cell(Lattice const &lattice) const;

private:
  Points points_;
};

// The images of `c` under `ops` that are distinct modulo the cell. Their count
// is the cluster's per-cell multiplicity, which is why no orbit in this module
// stores a multiplicity that could drift out of step with its images.
[[nodiscard]] std::vector<Cluster>
symmetry_images(Cluster const &c, Operations const &ops, double tol);

// The site-symmetry subgroup of `c`: the operations mapping it onto itself
// modulo the cell, each corrected by the translation that makes it permute the
// cluster's own points with zero net shift.
//
// Those corrected translations are not lattice-commensurate. SymmetryOperation
// permits that -- only the rotation is integral, and it is the input rotation
// unchanged -- but it does mean the result is a stabilizer, not a subgroup of
// the space group as stored elsewhere in this library.
[[nodiscard]] Operations site_symmetry(Cluster const &c, Operations const &ops,
                                       double tol);

// How `op` permutes the points of `c`: entry i is the index of the point that
// op maps point i onto. Only meaningful for an op from site_symmetry(c, ...),
// which is exactly where the configuration builder uses it -- comparing
// occupations through a precomputed permutation beats rebuilding and matching
// decorated geometry for every candidate assignment.
[[nodiscard]] std::vector<int> point_permutation(Cluster const &c,
                                                 SymmetryOperation const &op,
                                                 double tol);

// The first element of `reps` whose projected cluster is symmetry-equivalent to
// `c`, else the end iterator. Callers that only want a yes/no discard the
// iterator; the ones that update a match in place keep it.
template <std::ranges::forward_range R, class Proj = std::identity>
  requires std::convertible_to<
      std::invoke_result_t<Proj &, std::ranges::range_reference_t<R>>,
      Cluster const &>
[[nodiscard]] auto find_equivalent(R &&reps, Cluster const &c,
                                   Operations const &ops, double tol,
                                   Proj proj = {}) {
  return std::ranges::find_if(reps, [&](auto const &e) {
    Cluster const &rep = std::invoke(proj, e);
    return std::ranges::any_of(ops, [&](SymmetryOperation const &op) {
      return c.transformed(op).congruent(rep, tol);
    });
  });
}

} // namespace seitz::alloy

#pragma GCC visibility pop
