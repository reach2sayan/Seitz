#include "cif_corpus.hpp"
#include "helpers.hpp"
#include "oracle.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/data/spacegroup_symbols.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/io/cif.hpp>

#include <catch2/catch_test_macros.hpp>

#include <optional>
#include <ranges>
#include <string>
#include <vector>

using namespace seitz;
using seitz::test::errored;
using seitz::test::must;

namespace {

// The space group a `miscellaneous/cifs` filename names, for the files whose
// name is a symbol the tables carry. The rest are non-standard settings
// (`Bmab`, `Cmca`, `C-1`, `F1`) that no table lists; the count assertion below
// is what keeps that set from silently growing.
[[nodiscard]] std::optional<int> number_from_name(std::string const &name) {
  auto const hall = data::hall_from_hm_symbol<GroupFamily::space>(name);
  return hall ? std::optional{data::spacegroup_type(*hall).number}
              : std::nullopt;
}

[[nodiscard]] std::vector<io::CifStructure> structures_of(
    oracle::CifFile const &file) {
  return must(io::read_cif(file.text));
}

} // namespace

TEST_CASE("every per-space-group CIF reads", "[oracle][cif]") {
  auto const corpus = oracle::cif_corpus("miscellaneous");
  REQUIRE(corpus.size() == 209);

  int resolved = 0;
  for (auto const &file : corpus) {
    INFO(file.name);
    auto const structures = structures_of(file);
    REQUIRE(structures.size() == 1);
    auto const &structure = structures.front();
    CHECK(structure.cell.size() > 0);

    auto const expected = number_from_name(file.name);
    if (!expected) {
      continue;
    }
    ++resolved;
    // What the file itself says has to be the group its name promises. This
    // is a statement about the reader's symbol resolution, not about the
    // atoms: the arrangement often has a higher symmetry than the setting it
    // is described in, which the determination below is free to find.
    REQUIRE(structure.hall.has_value());
    CHECK(data::spacegroup_type(*structure.hall).number == *expected);
  }
  // The non-standard settings are a fixed, small set.
  CHECK(resolved >= 190);
}

TEST_CASE("the determination of a read cell matches the reference",
          "[oracle][cif]") {
  for (auto const &file : oracle::cif_corpus("miscellaneous")) {
    INFO(file.name);
    auto const structure = structures_of(file).front();
    auto const ours = test::dataset_of(structure.cell, io::kCifTolerance);
    REQUIRE(ours);
    auto const reference =
        oracle::reference_dataset(structure.cell, io::kCifTolerance.symprec);
    CHECK(data::spacegroup_type(ours->hall).number == reference.number);
  }
}

TEST_CASE("a written CIF reads back as the same structure", "[oracle][cif]") {
  for (auto const &file : oracle::cif_corpus("miscellaneous")) {
    INFO(file.name);
    auto const structure = structures_of(file).front();

    // P1: every atom written out, so the round trip is exact.
    auto const plain = must(io::read_cif(io::write_cif(structure.cell)));
    REQUIRE(plain.size() == 1);
    CHECK(plain.front().cell.size() == structure.cell.size());

    // Symmetrized: one representative per orbit plus the database operations,
    // which have to expand back to the same standardized cell.
    auto const analyzer = analysis::SymmetryAnalyzer::from_cell(
        structure.cell, io::kCifTolerance);
    auto const text = io::write_cif(analyzer);
    REQUIRE(text);
    auto const again = must(io::read_cif(*text));
    REQUIRE(again.size() == 1);
    REQUIRE(again.front().hall.has_value());
    CHECK(data::spacegroup_type(*again.front().hall).number ==
          data::spacegroup_type(must(analyzer.hall())).number);
    CHECK(again.front().cell.size() ==
          must(analyzer.standardized_cell()).size());
  }
}

TEST_CASE("the real-structure corpus reads", "[oracle][cif]") {
  auto const corpus = oracle::cif_corpus("database");
  REQUIRE(corpus.size() == 77);
  for (auto const &file : corpus) {
    INFO(file.name);
    auto const structures = structures_of(file);
    CHECK(!structures.empty());
    for (auto const &structure : structures) {
      CHECK(structure.cell.size() > 0);
      CHECK(structure.hall.has_value());
    }
  }
}
