// Broad oracle sweep (Phase 9c): run get_dataset over spglib's full reference
// corpus (one input cell per space group, ~230 cells) and diff the
// determination against spg_get_dataset — flushing latent bugs that the small
// hand-picked corpus misses. Also cross-checks against the expected space-group
// number embedded in each YAML.

#include "corpus.hpp"
#include "oracle.hpp"

#include <cppcrystal/dataset.hpp>

#include <catch2/catch_test_macros.hpp>

#include <string>

TEST_CASE("get_dataset matches spg_get_dataset across the reference corpus",
          "[oracle][corpus]") {
  double const symprec = 1e-5;
  auto const corpus = cppcrystal::oracle::load_corpus();
  INFO("corpus size " << corpus.size());
  REQUIRE(corpus.size() >= 200); // ~230 cells; guards against a bad data path

  std::size_t checked = 0;
  for (auto const &entry : corpus) {
    INFO("cell " << entry.name << " (expected SG "
                 << entry.space_group_number << ")");

    auto const got = cppcrystal::get_dataset(entry.cell, {symprec});
    REQUIRE(got);
    auto const ref = cppcrystal::oracle::reference_dataset(entry.cell, symprec);
    REQUIRE(ref.number != 0);

    CHECK(got->spacegroup_number == ref.number);
    CHECK(got->spacegroup_number == entry.space_group_number);
    CHECK(got->hall_number == ref.hall_number);
    CHECK(got->international_symbol == std::string_view(ref.international));
    CHECK(static_cast<int>(got->operations.size()) == ref.n_operations);
    CHECK(static_cast<int>(got->std_types.size()) == ref.n_std_atoms);
    ++checked;
  }
  CHECK(checked >= 200);
}
