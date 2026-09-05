#pragma once

#include <seitz/alloy/cluster.hpp>
#include <seitz/alloy/parent_lattice.hpp>
#include <seitz/alloy/site_basis.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/types.hpp>

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::alloy {

// A decorated subcluster: one basis function of the expansion. These index the
// correlation vector xi, hence the columns of every v-matrix. Multiplicity =
// |images|, and the v-matrix builder reads the same images.
struct ClusterFunction {
  Cluster cluster;
  std::vector<Cluster> images;

  [[nodiscard]] int multiplicity() const noexcept {
    return static_cast<int>(images.size());
  }
};

// One species assignment on a subcluster, with the size of its orbit under
// that subcluster's site symmetry. Occupations are species indices in [0, k),
// kept apart from the point-function indices in [0, k - 1).
struct Configuration {
  std::vector<int> occupation; // one species index per point of the subcluster
  int multiplicity = 0;
};

// One symmetry-distinct subcluster, carrying what the CVM entropy
//
//   S = -k_B sum_c k_c sum_j m_jc rho_jc ln rho_jc,   rho_.c = V_c . xi
//
// needs: the configurations j with multiplicities m_jc, the v-matrix V_c
// mapping xi onto cluster probabilities, and the Kikuchi-Barker coefficient
// k_c.
struct CvmCluster {
  Cluster sites; // undecorated: every point function is zero
  std::vector<Configuration> configurations;
  MatrixXd vmatrix; // rows = configurations, cols = ALL cluster functions
  double kikuchi_barker = 0.0;
  double diameter = 0.0;
};

// The Cluster Variation Method ingredients of a set of maximal clusters.
class Cvm {
public:
  // `maximal`: the maximal clusters in the parent's fractional frame with
  // species counts filled in, as ParentLattice::cluster_of produces. `basis`
  // must cover the parent's largest sublattice; it is read once, not kept.
  [[nodiscard]] static Result<Cvm> create(ParentLattice const &parent,
                                          std::span<Cluster const> maximal,
                                          SiteBasis const &basis,
                                          double tolerance = kClusterPrec);

  // The subclusters, ordered by (point count, diameter); the empty cluster is
  // first.
  [[nodiscard]] std::span<CvmCluster const> clusters() const noexcept {
    return clusters_;
  }
  // The cluster functions that index the correlation vector, in subcluster
  // order.
  [[nodiscard]] std::span<ClusterFunction const> functions() const noexcept {
    return functions_;
  }
  // Carried so a consumer can measure a subcluster or recover Cartesian
  // geometry without the parent alongside.
  [[nodiscard]] Lattice const &lattice() const noexcept { return lattice_; }

private:
  Cvm(Lattice lattice, std::vector<CvmCluster> clusters,
      std::vector<ClusterFunction> functions)
      : lattice_{std::move(lattice)}, clusters_{std::move(clusters)},
        functions_{std::move(functions)} {}

  Lattice lattice_;
  std::vector<CvmCluster> clusters_;
  std::vector<ClusterFunction> functions_;
};

} // namespace seitz::alloy

#pragma GCC visibility pop
