#pragma once

#include <cppcrystal/alloy/cluster.hpp>
#include <cppcrystal/alloy/parent_lattice.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/types.hpp>

#include <boost/container/flat_map.hpp>

#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace cppcrystal::alloy {

// One symmetry orbit of decorated clusters: a representative and its images
// distinct modulo the cell.
//
// The multiplicity is DERIVED from the images rather than stored beside them.
// It is their count by definition, so a stored copy would be one more thing to
// keep true, and callers compare an integer instead of a double.
struct Orbit {
  Cluster representative;
  std::vector<Cluster> images;
  double diameter = 0.0;

  [[nodiscard]] int multiplicity() const noexcept {
    return static_cast<int>(images.size());
  }
};

// The symmetry-distinct decorated clusters of a parent lattice, up to a
// per-body-order diameter cutoff.
//
// Immutable once generated, like every other container in this library: a pool
// is derived with subpool(), never edited in place.
class ClustersPool {
public:
  struct Options {
    // Largest diameter to keep, per body order. Keyed by order so a caller
    // writes {{2, 6.0}, {3, 4.5}} and the orders present are exactly the orders
    // enumerated -- no unused entry, no parallel arrays to keep in step, and no
    // way to name an order without a cutoff.
    boost::container::flat_map<int, double> radii;
    bool include_empty = true;  // the constant (order-0) cluster
    bool include_points = true; // the order-1 point clusters
    Anchors anchors = Anchors::inequivalent;
    double tolerance = kClusterPrec;
  };

  // Orbits come out ordered: the empty cluster, then the points, then each body
  // order ascending and, within an order, by increasing diameter.
  [[nodiscard]] static Result<ClustersPool> generate(ParentLattice const &parent,
                                                     Options const &options);

  [[nodiscard]] auto begin() const noexcept { return orbits_.begin(); }
  [[nodiscard]] auto end() const noexcept { return orbits_.end(); }
  [[nodiscard]] bool empty() const noexcept { return orbits_.empty(); }
  [[nodiscard]] Index size() const noexcept {
    return static_cast<Index>(orbits_.size());
  }
  [[nodiscard]] Orbit const &operator[](Index i) const noexcept {
    return orbits_[static_cast<std::size_t>(i)];
  }

  // The parent's lattice, carried along because every consumer of an orbit
  // needs it (diameters, Cartesian geometry) and would otherwise have to keep
  // the parent alive beside the pool.
  [[nodiscard]] Lattice const &lattice() const noexcept { return lattice_; }

  // A pool holding the orbits at `indices`, in the order given: what a
  // downstream consumer takes when it wants a subset rather than all of them.
  template <std::ranges::input_range R>
  [[nodiscard]] ClustersPool subpool(R &&indices) const {
    auto picked = std::forward<R>(indices) |
                  std::views::transform([this](auto i) -> Orbit const & {
                    return orbits_[static_cast<std::size_t>(i)];
                  });
    return ClustersPool{lattice_,
                        std::vector<Orbit>(picked.begin(), picked.end())};
  }

private:
  ClustersPool(Lattice lattice, std::vector<Orbit> orbits)
      : lattice_{std::move(lattice)}, orbits_{std::move(orbits)} {}

  Lattice lattice_;
  std::vector<Orbit> orbits_;
};

} // namespace cppcrystal::alloy

#pragma GCC visibility pop
