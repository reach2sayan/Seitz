// Alloy cluster enumeration and the Cluster Variation Method.
//
// The reference values are ATAT's, on the binary fcc lattice: the pair and
// triangle multiplicities, and every ingredient the CVM builds on the
// nearest-neighbour tetrahedron. They are textbook numbers for fcc, so they
// pin the port without needing an external oracle.
#include <cppcrystal/alloy/clusters_pool.hpp>
#include <cppcrystal/alloy/cvm.hpp>
#include <cppcrystal/alloy/parent_lattice.hpp>
#include <cppcrystal/alloy/site_basis.hpp>

#include "helpers.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <ranges>
#include <vector>

using namespace cppcrystal;
using namespace cppcrystal::alloy;
using cppcrystal::test::errored;
using cppcrystal::test::must;
using Catch::Approx;

namespace {

// Aluminium and titanium stand in for ATAT's "A,B"; only the count matters to
// the mathematics, but the parent lattice speaks atomic numbers like the rest
// of the library.
constexpr int kA = 13;
constexpr int kB = 22;
constexpr int kC = 24;

// The fcc primitive cell with a = 1: nearest neighbours at 1/sqrt(2), second
// neighbours at 1. The matrix is symmetric, so the rows-versus-columns
// question the lattice file format leaves open does not arise here.
[[nodiscard]] Lattice fcc_primitive() {
  Matrix3d basis;
  basis << 0.0, 0.5, 0.5, //
      0.5, 0.0, 0.5,      //
      0.5, 0.5, 0.0;
  return *Lattice::from_basis(basis);
}

[[nodiscard]] ParentLattice binary_fcc(Species allowed = Species{kA, kB}) {
  std::array const sites{SiteSpec{Vector3d::Zero(), std::move(allowed)}};
  return must(ParentLattice::from_sites(fcc_primitive(), sites));
}

// The nearest-neighbour tetrahedron, in Cartesian coordinates: the origin and
// three of its twelve nearest neighbours.
[[nodiscard]] Cluster nn_tetrahedron(ParentLattice const &parent) {
  std::array<Vector3d, 4> const corners{
      Vector3d{0.0, 0.0, 0.0}, Vector3d{0.0, 0.5, 0.5},
      Vector3d{0.5, 0.0, 0.5}, Vector3d{0.5, 0.5, 0.0}};
  return must(parent.cluster_of(corners));
}

[[nodiscard]] Cvm tetrahedron_cvm(ParentLattice const &parent) {
  std::array const maximal{nn_tetrahedron(parent)};
  return must(Cvm::create(parent, maximal,
                          SiteBasis::trigonometric(parent.max_species())));
}

// The orbits of a given body order.
[[nodiscard]] std::vector<Orbit const *> of_order(ClustersPool const &pool,
                                                  int order) {
  std::vector<Orbit const *> found;
  for (Orbit const &orbit : pool) {
    if (orbit.representative.size() == order) {
      found.push_back(&orbit);
    }
  }
  return found;
}

} // namespace

TEST_CASE("alloy: the trigonometric site basis is the Ising spin on a binary",
          "[alloy]") {
  auto const basis = SiteBasis::trigonometric(2);
  CHECK(basis.max_species() == 2);
  CHECK(basis[2, 0, 0] == Approx(-1.0));
  CHECK(basis[2, 0, 1] == Approx(1.0));
  // The v-matrix's per-site normalization: theta_0(0)^2 + theta_0(1)^2.
  CHECK(basis.sum_of_squares(2, 0) == Approx(2.0));
}

TEST_CASE("alloy: fcc pair orbits carry the textbook multiplicities",
          "[alloy]") {
  auto const parent = binary_fcc();
  REQUIRE(parent.operations().size() == 48);

  ClustersPool::Options options;
  options.radii = {{2, 1.05}}; // reaches both the 1st and 2nd neighbour shells
  auto const pool = must(ClustersPool::generate(parent, options));

  REQUIRE(of_order(pool, 0).size() == 1); // the empty cluster
  REQUIRE(of_order(pool, 1).size() == 1); // the point
  auto const pairs = of_order(pool, 2);
  REQUIRE(pairs.size() == 2);

  // Ordered by increasing diameter within the body order.
  CHECK(pairs[0]->diameter == Approx(0.70710678).epsilon(1e-6));
  CHECK(pairs[0]->multiplicity() == 6);
  CHECK(pairs[1]->diameter == Approx(1.0).epsilon(1e-6));
  CHECK(pairs[1]->multiplicity() == 3);
}

TEST_CASE("alloy: the fcc nearest-neighbour triangle has multiplicity 8",
          "[alloy]") {
  auto const parent = binary_fcc();
  ClustersPool::Options options;
  options.radii = {{2, 0.8}, {3, 0.8}}; // nearest-neighbour geometry only
  auto const pool = must(ClustersPool::generate(parent, options));

  REQUIRE(of_order(pool, 2).size() == 1);
  auto const triangles = of_order(pool, 3);
  REQUIRE(triangles.size() == 1);
  CHECK(triangles[0]->multiplicity() == 8);
  CHECK(triangles[0]->diameter == Approx(0.70710678).epsilon(1e-6));
}

TEST_CASE("alloy: with no radii the pool is the empty and point clusters",
          "[alloy]") {
  auto const parent = binary_fcc();
  auto const pool = must(ClustersPool::generate(parent, {}));

  REQUIRE(pool.size() == 2);
  CHECK(pool[0].representative.empty());
  CHECK(pool[0].multiplicity() == 1);
  CHECK(pool[1].representative.size() == 1);
  CHECK(pool[1].multiplicity() == 1);
  CHECK(pool[1].representative[0].function == 0);
  CHECK(pool[1].representative[0].species == 2);
}

TEST_CASE("alloy: seeding on inequivalent sites gives the same orbits as "
          "seeding on all of them",
          "[alloy]") {
  // The fcc lattice as its 4-atom conventional cell, where the active sites
  // are one crystallographic orbit rather than one site. The 1-atom primitive
  // fixture above cannot tell the two anchor strategies apart.
  std::array const sites{
      SiteSpec{Vector3d{0.0, 0.0, 0.0}, Species{kA, kB}},
      SiteSpec{Vector3d{0.0, 0.5, 0.5}, Species{kA, kB}},
      SiteSpec{Vector3d{0.5, 0.0, 0.5}, Species{kA, kB}},
      SiteSpec{Vector3d{0.5, 0.5, 0.0}, Species{kA, kB}}};
  auto const parent =
      must(ParentLattice::from_sites(Lattice{Matrix3d::Identity()}, sites));

  REQUIRE(parent.anchors(Anchors::all).size() == 4);
  REQUIRE(parent.anchors(Anchors::inequivalent).size() == 1);

  ClustersPool::Options options;
  options.radii = {{2, 1.05}, {3, 0.8}};

  options.anchors = Anchors::inequivalent;
  auto const fast = must(ClustersPool::generate(parent, options));
  options.anchors = Anchors::all;
  auto const exhaustive = must(ClustersPool::generate(parent, options));

  REQUIRE(fast.size() == exhaustive.size());
  for (auto const &[quick, full] : std::views::zip(fast, exhaustive)) {
    CHECK(quick.representative.size() == full.representative.size());
    CHECK(quick.multiplicity() == full.multiplicity());
    CHECK(quick.diameter == Approx(full.diameter).margin(1e-9));
  }
}

TEST_CASE("alloy: fcc tetrahedron CVM subclusters and multiplicities",
          "[alloy]") {
  auto const parent = binary_fcc();
  auto const cvm = tetrahedron_cvm(parent);

  // Empty, point, pair, triplet, quadruplet; a binary decorates each exactly
  // one way, so there are as many cluster functions as subclusters.
  REQUIRE(cvm.clusters().size() == 5);
  REQUIRE(cvm.functions().size() == 5);

  std::vector<int> sizes;
  std::vector<int> multiplicities;
  for (ClusterFunction const &function : cvm.functions()) {
    sizes.push_back(function.cluster.size());
    multiplicities.push_back(function.multiplicity());
  }
  CHECK(sizes == std::vector<int>{0, 1, 2, 3, 4});
  CHECK(multiplicities == std::vector<int>{1, 1, 6, 8, 2});
}

TEST_CASE("alloy: fcc tetrahedron configurations are binomial", "[alloy]") {
  auto const parent = binary_fcc();
  auto const cvm = tetrahedron_cvm(parent);

  std::vector<std::size_t> counts;
  std::vector<std::vector<int>> weights;
  for (CvmCluster const &cluster : cvm.clusters()) {
    counts.push_back(cluster.configurations.size());
    std::vector<int> here;
    for (Configuration const &configuration : cluster.configurations) {
      here.push_back(configuration.multiplicity);
    }
    weights.push_back(std::move(here));
  }
  CHECK(counts == std::vector<std::size_t>{1, 2, 3, 4, 5});
  CHECK(weights[2] == std::vector<int>{1, 2, 1});
  CHECK(weights[3] == std::vector<int>{1, 3, 3, 1});
  CHECK(weights[4] == std::vector<int>{1, 4, 6, 4, 1});
}

TEST_CASE("alloy: fcc Kikuchi-Barker coefficients and the pair v-matrix",
          "[alloy]") {
  auto const parent = binary_fcc();
  auto const cvm = tetrahedron_cvm(parent);

  std::vector<double> coefficients;
  for (CvmCluster const &cluster : cvm.clusters()) {
    coefficients.push_back(cluster.kikuchi_barker);
  }
  CHECK(coefficients == std::vector<double>{0, 5, -1, 0, 1});

  // Rows are the pair's three configurations, columns all five cluster
  // functions; only the empty, point and pair columns can overlap a pair.
  MatrixXd const &v = cvm.clusters()[2].vmatrix;
  REQUIRE(v.rows() == 3);
  REQUIRE(v.cols() == 5);
  constexpr double expected[3][5] = {{0.25, -0.5, 0.25, 0, 0},
                                     {0.25, 0.0, -0.25, 0, 0},
                                     {0.25, 0.5, 0.25, 0, 0}};
  for (Index row = 0; row < 3; ++row) {
    for (Index column = 0; column < 5; ++column) {
      CHECK(v(row, column) ==
            Approx(expected[row][column]).margin(1e-9));
    }
  }

  // The empty subcluster's probability is identically one: its single row
  // picks out the empty correlation and nothing else.
  MatrixXd const &constant = cvm.clusters()[0].vmatrix;
  REQUIRE(constant.rows() == 1);
  CHECK(constant(0, 0) == Approx(1.0).margin(1e-9));
  CHECK(constant.rightCols(4).cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-9));
}

TEST_CASE("alloy: a ternary point cluster carries two cluster functions",
          "[alloy]") {
  auto const parent = binary_fcc(Species{kA, kB, kC});
  auto const cvm = tetrahedron_cvm(parent);

  int point_functions = 0;
  for (ClusterFunction const &function : cvm.functions()) {
    if (function.cluster.size() == 1) {
      ++point_functions;
    }
  }
  CHECK(point_functions == 2);
  CHECK(cvm.clusters()[1].configurations.size() == 3); // A, B or C on the point
}

TEST_CASE("alloy: a parent lattice validates its sublattices", "[alloy]") {
  Cell const cell{fcc_primitive(), Positions::Zero(1, 3), Types{1}};
  // Sublattice id 1 with only one entry in the species table.
  CHECK(errored([&] {
    return ParentLattice::create(cell, std::vector<Species>{Species{kA, kB}});
  }));
  CHECK(errored([&] {
    return ParentLattice::create(Cell{fcc_primitive(), Positions::Zero(1, 3),
                                      Types{0}},
                                 std::vector<Species>{Species{kB, kA}});
  }));
  CHECK(errored([&] {
    return ParentLattice::from_sites(fcc_primitive(), std::span<SiteSpec>{});
  }));
}

TEST_CASE("alloy: a spectator sublattice carries no cluster point", "[alloy]") {
  std::array const sites{SiteSpec{Vector3d::Zero(), Species{kA, kB}},
                         SiteSpec{Vector3d{0.25, 0.25, 0.25}, Species{kC}}};
  auto const parent = must(ParentLattice::from_sites(fcc_primitive(), sites));

  CHECK(parent.max_species() == 2);
  auto const anchors = parent.anchors(Anchors::all);
  REQUIRE(anchors.size() == 1);
  CHECK(anchors[0].species == 2);

  auto const pool = must(ClustersPool::generate(parent, {}));
  CHECK(pool.size() == 2);
}

TEST_CASE("alloy: from_sites derives the sublattice ids create expects",
          "[alloy]") {
  // The species lists are written in different orders and one is repeated:
  // sorting makes "same chemistry" the thing that decides, and the repeat
  // collapses onto one id.
  std::array const sites{SiteSpec{Vector3d::Zero(), Species{kB, kA}},
                         SiteSpec{Vector3d{0.25, 0.25, 0.25}, Species{kA, kB}},
                         SiteSpec{Vector3d{0.5, 0.5, 0.5}, Species{kC}}};
  auto const parent = must(ParentLattice::from_sites(fcc_primitive(), sites));

  CHECK(parent.cell().types() == Types{0, 0, 1});
  CHECK(parent.species_at(0) == 2);
  CHECK(parent.species_at(2) == 1);
}

TEST_CASE("alloy: a cluster round-trips through Cell and subpool", "[alloy]") {
  auto const parent = binary_fcc();
  ClustersPool::Options options;
  options.radii = {{2, 0.8}};
  auto const pool = must(ClustersPool::generate(parent, options));
  REQUIRE(pool.size() == 3); // empty, point, nearest-neighbour pair

  Cell const pair = pool[2].representative.as_cell(pool.lattice());
  CHECK(pair.size() == 2);
  CHECK(pair.position(0).isApprox(pool[2].representative[0].position));
  CHECK(pair.types() == Types{0, 0}); // both points carry function 0

  // subpool keeps the order it is given, not the pool's.
  auto const picked = pool.subpool(std::array{2, 0});
  REQUIRE(picked.size() == 2);
  CHECK(picked[0].representative.size() == 2);
  CHECK(picked[1].representative.empty());
}
