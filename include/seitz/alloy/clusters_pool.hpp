#pragma once

#include <seitz/alloy/cluster.hpp>
#include <seitz/alloy/parent_lattice.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/types.hpp>

#include <boost/container/flat_map.hpp>

#include <cstddef>
#include <ranges>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::alloy {

// One symmetry orbit of decorated clusters: a representative and its images
// distinct modulo the cell. Multiplicity = |images|, derived not stored.
struct Orbit {
  Cluster representative;
  std::vector<Cluster> images;
  double diameter = 0.0;

  [[nodiscard]] int multiplicity() const noexcept {
    return static_cast<int>(images.size());
  }
};

// The symmetry-distinct decorated clusters of a parent lattice, up to a
// per-body-order diameter cutoff. Immutable: derive with subpool().
class ClustersPool {
public:
  struct Options {
    // Largest diameter kept, per body order: {{2, 6.0}, {3, 4.5}}. The keys
    // are exactly the orders enumerated.
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

  // The parent's lattice: every consumer of an orbit needs it (diameters,
  // Cartesian geometry) without holding the parent alive.
  [[nodiscard]] Lattice const &lattice() const noexcept { return lattice_; }

  // A pool holding the orbits at `indices`, in the order given.
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

} // namespace seitz::alloy

#pragma GCC visibility pop
