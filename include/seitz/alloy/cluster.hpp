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

// Coincidence tolerance for cluster points, fractional. Geometry is exact
// symmetry acting on exact lattice positions, so the only residual is the 1/n
// averaged translation of site_symmetry (~1e-16); the value is ATAT's global
// tolerance, which reproduces its multiplicities.
inline constexpr double kClusterPrec = 1e-3;

// One cluster point: position, the species count k its sublattice admits, and
// the point function decorating it. `function` indexes the site basis, in
// [0, k - 1); a species assignment lives in Configuration::occupation, not
// here -- the two ranges differ and sharing a field conflates them.
struct ClusterPoint {
  Vector3d position{Vector3d::Zero()}; // fractional, in the parent's frame
  int species = 2;
  int function = 0;
};

// A decorated cluster: a set of points, each carrying its point function.
// small_vector<_, 4>: the nearest-neighbour tetrahedron is the CVM's workhorse
// maximal cluster and enumeration builds millions of these.
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

  // max_{i,j} |r_i - r_j| in Cartesian: the cluster length compared by the
  // radius cutoffs and the subcluster ordering.
  [[nodiscard]] double diameter(Lattice const &lattice) const;

  // Equal modulo a lattice translation: some whole-cluster shift maps every
  // point of *this onto a point of `other` with the same function. Fractional
  // frame, so "modulo the cell" is "modulo 1" -- a rounding test, not a metric
  // one. `species` is deliberately not compared (as in ATAT): the function
  // index identifies the basis element, the sublattice size does not.
  [[nodiscard]] bool congruent(Cluster const &other, double tol) const;

  // The cluster as a Cell. Atom types are the point FUNCTIONS, so a symmetry
  // determination on it sees the decorated cluster, not the bare geometry.
  [[nodiscard]] Cell as_cell(Lattice const &lattice) const;

private:
  Points points_;
};

// Images of `c` under `ops` distinct modulo the cell; their count is the
// per-cell multiplicity, hence never stored beside them.
[[nodiscard]] std::vector<Cluster>
symmetry_images(Cluster const &c, Operations const &ops, double tol);

// Stabilizer of `c`: the ops mapping it onto itself modulo the cell, each
// corrected by the translation that permutes c's own points with zero net
// shift. Those translations are not lattice-commensurate, so this is a
// stabilizer, not a subgroup of the space group.
[[nodiscard]] Operations site_symmetry(Cluster const &c, Operations const &ops,
                                       double tol);

// The permutation `op` induces on the points of `c`: entry i is the image of
// point i. Meaningful only for op in site_symmetry(c, ...); the configuration
// builder compares occupations through it.
[[nodiscard]] std::vector<int> point_permutation(Cluster const &c,
                                                 SymmetryOperation const &op,
                                                 double tol);

// First element of `reps` whose projected cluster is symmetry-equivalent to
// `c`, else the end iterator.
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
