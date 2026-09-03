// Robustness sweep (Phase 9c): the space-group determination must be stable
// under sub-symprec position/lattice jitter and must return the same space group
// for an integer supercell. Exercises the tolerance-retry loops over a sample of
// the reference corpus.

#include "corpus.hpp"

#include <cppcrystal/analysis/symmetry_analyzer.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

#include <random>
#include <vector>

namespace {

using namespace cppcrystal;

// Perturb every atom by a random Cartesian displacement bounded well below
// symprec (so symmetry-related atoms still coincide within symprec).
Cell jitter(Cell const &cell, std::mt19937 &rng, double cart_amplitude) {
  Matrix3d const lattice_inv = cell.lattice().matrix().inverse();
  std::uniform_real_distribution<double> dist(-cart_amplitude, cart_amplitude);
  Positions positions = cell.positions();
  for (Index i = 0; i < cell.size(); ++i) {
    Vector3d const cart(dist(rng), dist(rng), dist(rng));
    positions.row(i) += (lattice_inv * cart).transpose();
  }
  return Cell(Lattice{cell.lattice().matrix()}, positions, cell.types());
}

// A diagonal n x n x n supercell (same space group as the input).
Cell supercell(Cell const &cell, int n) {
  auto const nd = static_cast<double>(n);
  Matrix3d const lattice = cell.lattice().matrix() * nd;
  Index const atoms = cell.size();
  Index const reps = static_cast<Index>(n) * n * n;
  Positions positions(atoms * reps, 3);
  std::vector<int> types(static_cast<std::size_t>(atoms * reps));
  Index out = 0;
  for (Index a = 0; a < atoms; ++a) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          Vector3d const frac =
              (cell.position(a) + Vector3d(i, j, k)) / nd;
          positions.row(out) = frac.transpose();
          types[static_cast<std::size_t>(out)] = cell.type(a);
          ++out;
        }
      }
    }
  }
  return Cell(Lattice{lattice}, positions, types);
}

} // namespace

TEST_CASE("get_dataset is stable under sub-symprec jitter",
          "[oracle][robustness]") {
  double const symprec = 1e-5;
  std::mt19937 rng(0xC0FFEE);
  auto const corpus = cppcrystal::oracle::load_corpus();
  REQUIRE(corpus.size() >= 200);

  for (std::size_t i = 0; i < corpus.size(); i += 13) {
    auto const &entry = corpus[i];
    INFO("cell " << entry.name << " (SG " << entry.space_group_number << ")");
    Cell const jittered = jitter(entry.cell, rng, 0.15 * symprec);
    auto const got = cppcrystal::test::dataset_of(jittered, {symprec});
    REQUIRE(got);
    CHECK(data::spacegroup_type(got->hall).number == entry.space_group_number);
  }
}

TEST_CASE("get_dataset returns the same space group for a 2x2x2 supercell",
          "[oracle][robustness]") {
  double const symprec = 1e-5;
  auto const corpus = cppcrystal::oracle::load_corpus();
  REQUIRE(corpus.size() >= 200);

  for (std::size_t i = 0; i < corpus.size(); i += 17) {
    auto const &entry = corpus[i];
    INFO("cell " << entry.name << " (SG " << entry.space_group_number << ")");
    Cell const sc = supercell(entry.cell, 2);
    auto const got = cppcrystal::test::dataset_of(sc, {symprec});
    REQUIRE(got);
    CHECK(data::spacegroup_type(got->hall).number == entry.space_group_number);
  }
}
