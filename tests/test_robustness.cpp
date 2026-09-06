// Robustness sweep (Phase 9c): the space-group determination must be stable
// under sub-symprec position/lattice jitter and must return the same space
// group for an integer supercell. Exercises the tolerance-retry loops over a
// sample of the reference corpus.

#include "corpus.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

#include <random>
#include <vector>

using namespace seitz;

TEST_CASE("get_dataset is stable under sub-symprec jitter",
          "[oracle][robustness]") {
  double const symprec = 1e-5;
  std::mt19937 rng(0xC0FFEE);
  auto const corpus = seitz::oracle::load_corpus();
  REQUIRE(corpus.size() >= 200);

  for (std::size_t i = 0; i < corpus.size(); i += 13) {
    auto const &entry = corpus[i];
    INFO("cell " << entry.name << " (SG " << entry.space_group_number << ")");
    // Well below symprec, so symmetry-related atoms still coincide.
    Cell const jittered = test::jitter(entry.cell, rng, 0.15 * symprec);
    auto const got = seitz::test::dataset_of(jittered, {symprec});
    REQUIRE(got);
    CHECK(data::spacegroup_type(got->hall).number == entry.space_group_number);
  }
}

TEST_CASE("get_dataset returns the same space group for a 2x2x2 supercell",
          "[oracle][robustness]") {
  double const symprec = 1e-5;
  auto const corpus = seitz::oracle::load_corpus();
  REQUIRE(corpus.size() >= 200);

  for (std::size_t i = 0; i < corpus.size(); i += 17) {
    auto const &entry = corpus[i];
    INFO("cell " << entry.name << " (SG " << entry.space_group_number << ")");
    Cell const sc = test::must(entry.cell.supercell(Vector3i::Constant(2)));
    auto const got = seitz::test::dataset_of(sc, {symprec});
    REQUIRE(got);
    CHECK(data::spacegroup_type(got->hall).number == entry.space_group_number);
  }
}
