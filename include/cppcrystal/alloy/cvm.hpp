#pragma once

#include <cppcrystal/alloy/cluster.hpp>
#include <cppcrystal/alloy/parent_lattice.hpp>
#include <cppcrystal/alloy/site_basis.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/types.hpp>

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace cppcrystal::alloy {

// A decorated subcluster: one basis function of the expansion. These index the
// correlation vector, and so the columns of every v-matrix.
//
// Like Orbit, it stores its images and derives the multiplicity from them --
// which here also removes the reference implementation's repeated
// symmetry_images() call inside the v-matrix builder, since the images the
// count comes from are the ones that builder needs.
struct ClusterFunction {
  Cluster cluster;
  std::vector<Cluster> images;

  [[nodiscard]] int multiplicity() const noexcept {
    return static_cast<int>(images.size());
  }
};

// One species assignment on a subcluster, with the size of its orbit under that
// subcluster's own site symmetry.
//
// The occupation is its own field rather than being stuffed into the cluster's
// point functions, as the reference implementation does. Those are different
// quantities over different ranges -- a basis function index in [0, k-1) and a
// species index in [0, k) -- and sharing one field is what forced its v-matrix
// builder to read two different meanings off the same name in one loop body.
struct Configuration {
  std::vector<int> occupation; // one species index per point of the subcluster
  int multiplicity = 0;
};

// One symmetry-distinct subcluster, carrying everything the CVM entropy
//
//   S = -k_B sum_c k_c sum_j m_jc rho_jc ln rho_jc,   rho_.c = V_c . xi
//
// needs from it: its configurations j with their multiplicities m_jc, the
// v-matrix V_c that turns the correlation vector xi into cluster probabilities,
// and its Kikuchi-Barker coefficient k_c.
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
  // `maximal` are the maximal clusters in the parent's fractional frame with
  // each point's species count filled in -- which is exactly what
  // ParentLattice::cluster_of produces from Cartesian points.
  //
  // `basis` must cover the parent's largest sublattice; it is taken by value
  // because the v-matrices are built from it once and it is not kept.
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
  // geometry without holding the parent lattice alongside.
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

} // namespace cppcrystal::alloy

#pragma GCC visibility pop
