#include "helpers.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/data/element_data.hpp>
#include <seitz/data/spacegroup_symbols.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/io/cif.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <format>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

using namespace seitz;
using seitz::test::errored;
using seitz::test::must;

namespace {

[[nodiscard]] io::CifBlock one_block(std::string_view text) {
  auto blocks = must(io::parse_cif(text));
  REQUIRE(blocks.size() == 1);
  return std::move(blocks.front());
}

// A whole column as a vector, so it compares in one CHECK.
[[nodiscard]] std::vector<std::string> rows(io::CifBlock const &block,
                                            std::string_view tag) {
  auto const column = block.column(tag);
  REQUIRE(column.has_value());
  return {column->begin(), column->end()};
}

} // namespace

TEST_CASE("the block name and scalar items come back as one-row columns",
          "[cif]") {
  auto const block = one_block("data_nacl\n"
                               "_cell_length_a   5.64\n"
                               "_chemical_name   'sodium chloride'\n");
  CHECK(block.name == "nacl");
  CHECK(block.value("_cell_length_a") == "5.64");
  CHECK(block.value("_chemical_name") == "sodium chloride");
  CHECK(rows(block, "_cell_length_a") == std::vector<std::string>{"5.64"});
  CHECK(block.loops.empty());
}

TEST_CASE("tags are matched case-insensitively", "[cif]") {
  auto const block = one_block("data_x\n_Cell_Length_A 5.64\n");
  CHECK(block.value("_cell_length_a") == "5.64");
  CHECK(!block.value("_cell_length_b").has_value());
}

TEST_CASE("comments are skipped wherever they start", "[cif]") {
  auto const block = one_block("# a whole-line comment\n"
                               "data_x\n"
                               "_a 1  # trailing\n"
                               "_b 2\n");
  CHECK(block.value("_a") == "1");
  CHECK(block.value("_b") == "2");
}

TEST_CASE("quotes hold the delimiter when it is not followed by space",
          "[cif]") {
  auto const block = one_block("data_x\n"
                               "_a 'it's here'\n"
                               "_b \"say \"\"what\"\n"
                               "_c 'trailing space '\n");
  CHECK(block.value("_a") == "it's here");
  CHECK(block.value("_b") == "say \"\"what");
  CHECK(block.value("_c") == "trailing space ");
}

TEST_CASE("a semicolon text field spans lines and keeps inner semicolons",
          "[cif]") {
  auto const block = one_block("data_x\n"
                               "_note\n"
                               ";first line; still going\n"
                               "second line\n"
                               ";\n"
                               "_after 7\n");
  auto const note = block.value("_note");
  REQUIRE(note.has_value());
  CHECK(note->starts_with("first line; still going"));
  CHECK(note->find("second line") != std::string_view::npos);
  CHECK(block.value("_after") == "7");
}

TEST_CASE("the two null tokens read as an absent value", "[cif]") {
  auto const block = one_block("data_x\n_a ?\n_b .\n_c 0\n");
  CHECK(!block.value("_a").has_value());
  CHECK(!block.value("_b").has_value());
  CHECK(block.value("_c") == "0");
  // Absent as a value, present as a column: the file did say something.
  CHECK(block.column("_a").has_value());
}

TEST_CASE("a loop becomes one column per tag", "[cif]") {
  auto const block = one_block("data_x\n"
                               "loop_\n"
                               "_atom_site_label\n"
                               "_atom_site_fract_x\n"
                               "Na1 0.0\n"
                               "Cl1 0.5\n"
                               "_after 7\n");
  CHECK(rows(block, "_atom_site_label") ==
        std::vector<std::string>{"Na1", "Cl1"});
  CHECK(rows(block, "_atom_site_fract_x") ==
        std::vector<std::string>{"0.0", "0.5"});
  CHECK(block.value("_after") == "7");
  REQUIRE(block.loops.size() == 1);
  CHECK(block.loops.front() == std::vector<std::string>{"_atom_site_label",
                                                        "_atom_site_fract_x"});
}

TEST_CASE("two loops in one block stay separate", "[cif]") {
  auto const block = one_block("data_x\n"
                               "loop_ _symop 'x,y,z' '-x,-y,-z'\n"
                               "LOOP_ _label _sym Na1 Na\n");
  CHECK(rows(block, "_symop") == std::vector<std::string>{"x,y,z", "-x,-y,-z"});
  CHECK(rows(block, "_label") == std::vector<std::string>{"Na1"});
  CHECK(block.loops.size() == 2);
}

TEST_CASE("a loop mixes quoted and unquoted values", "[cif]") {
  auto const block = one_block("data_x\n"
                               "loop_ _a _b\n"
                               "1 'one'\n"
                               "2 \"two\"\n");
  CHECK(rows(block, "_b") == std::vector<std::string>{"one", "two"});
}

TEST_CASE("a document carries every data block in file order", "[cif]") {
  auto const blocks = must(io::parse_cif("data_first\n_a 1\n"
                                         "data_second\n_a 2\n"));
  REQUIRE(blocks.size() == 2);
  CHECK(blocks[0].name == "first");
  CHECK(blocks[1].name == "second");
  CHECK(blocks[1].value("_a") == "2");
}

TEST_CASE("a ragged loop reports the line of its own loop_ keyword", "[cif]") {
  auto const text = std::string_view{"data_x\n"
                                     "_a 1\n"
                                     "loop_ _p _q\n"
                                     "1 2 3\n"};
  CHECK(errored([&] { return io::parse_cif(text); }));
  auto const position = leaf::try_handle_all(
      [&]() -> Result<std::pair<std::int64_t, std::int64_t>> {
        BOOST_LEAF_AUTO(blocks, io::parse_cif(text));
        (void)blocks;
        return std::pair<std::int64_t, std::int64_t>{0, 0};
      },
      [](e_cif_syntax const &e) {
        return std::pair{e.line, e.column};
      },
      [](leaf::error_info const &) {
        return std::pair<std::int64_t, std::int64_t>{-1, -1};
      });
  CHECK(position == std::pair<std::int64_t, std::int64_t>{3, 1});
}

TEST_CASE("text that stops making sense reports where", "[cif]") {
  // The value never closes, so the item never completes.
  CHECK(errored([] { return io::parse_cif("data_x\n_a 'unterminated\n"); }));
  CHECK(errored([] { return io::parse_cif("not a cif at all\n"); }));
}

TEST_CASE("an empty or comment-only document has no blocks", "[cif]") {
  CHECK(must(io::parse_cif("")).empty());
  CHECK(must(io::parse_cif("   \n# nothing here\n")).empty());
}

// ---- the reader ------------------------------------------------------------

namespace {

// The six cell parameters every reader test needs, written once.
[[nodiscard]] std::string cubic_block(std::string_view body, double edge = 4.0,
                                      std::string_view name = "x") {
  return std::format("data_{0}\n"
                     "_cell_length_a {1}\n_cell_length_b {1}\n"
                     "_cell_length_c {1}\n"
                     "_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                     "_cell_angle_gamma 90\n{2}",
                     name, edge, body);
}

// Rocksalt as a file writes it; each test varies the symmetry part.
[[nodiscard]] std::string nacl_block(std::string_view symmetry) {
  return cubic_block(std::string{symmetry} +
                         "loop_\n"
                         "_atom_site_label _atom_site_type_symbol\n"
                         "_atom_site_fract_x _atom_site_fract_y "
                         "_atom_site_fract_z\n"
                         "Na1 Na 0.0 0.0 0.0\n"
                         "Cl1 Cl 0.5 0.5 0.5\n",
                     5.64, "nacl");
}

// The face-centring translations: enough to build the conventional cell.
constexpr std::string_view kFccSymops =
    "loop_ _space_group_symop_operation_xyz\n"
    "'x,y,z' 'x,y+1/2,z+1/2' 'x+1/2,y,z+1/2' 'x+1/2,y+1/2,z'\n";

[[nodiscard]] io::CifStructure one_structure(std::string_view text) {
  auto structures = must(io::read_cif(text));
  REQUIRE(structures.size() == 1);
  return std::move(structures.front());
}

} // namespace

TEST_CASE("a symop loop expands the asymmetric unit", "[cif]") {
  auto const structure = one_structure(nacl_block(kFccSymops));
  CHECK(structure.name == "nacl");
  CHECK(structure.cell.size() == 8);
  CHECK(structure.labels.size() == 8);
  CHECK(structure.collapsed.empty());
  // Every atom of a face-centred pair is the same species as its parent.
  CHECK(std::ranges::count(structure.cell.types(),
                           *data::atomic_number("Na")) == 4);
}

TEST_CASE("a named space group expands the same cell as its operations",
          "[cif]") {
  auto const from_symops = one_structure(nacl_block(kFccSymops));
  for (auto const naming : std::array<std::string_view, 4>{
           "_space_group_name_H-M_alt 'F m -3 m'\n",
           "_symmetry_space_group_name_H-M 'Fm-3m'\n",
           "_space_group_IT_number 225\n",
           "_space_group_name_Hall '-F 4 2 3'\n"}) {
    INFO(naming);
    auto const structure = one_structure(nacl_block(naming));
    REQUIRE(structure.hall.has_value());
    CHECK(data::spacegroup_type(*structure.hall).number == 225);
    // 48 operations x 2 sites, deduplicated back to the 8 distinct atoms.
    CHECK(structure.cell.size() == from_symops.cell.size());
  }
}

TEST_CASE("a block with no symmetry at all is read as P1", "[cif]") {
  auto const structure = one_structure(nacl_block(""));
  CHECK(!structure.hall.has_value());
  CHECK(structure.cell.size() == 2);
}

TEST_CASE("a shared site collapses to its majority species", "[cif]") {
  auto const text =
      cubic_block(
                  "loop_\n"
                  "_atom_site_label _atom_site_type_symbol\n"
                  "_atom_site_fract_x _atom_site_fract_y _atom_site_fract_z\n"
                  "_atom_site_occupancy\n"
                  "Na1 Na 0.0 0.0 0.0 0.6\n"
                  "K1  K  0.0 0.0 0.0 0.4\n"
                  "Cl1 Cl 0.5 0.5 0.5 1.0\n");
  auto const structure = one_structure(text);
  CHECK(structure.cell.size() == 2);
  CHECK(structure.cell.types() ==
        Types{*data::atomic_number("Na"), *data::atomic_number("Cl")});
  REQUIRE(structure.collapsed.size() == 1);
  CHECK(structure.collapsed.front().kept == "Na1");
  CHECK(structure.collapsed.front().dropped == std::vector<std::string>{"K1"});
}

TEST_CASE("an over-specified asymmetric unit reads as one orbit", "[cif]") {
  // The second row is the first's image under the stated operation, which a
  // redundantly written asymmetric unit does list.
  auto const text =
      cubic_block(
                  "loop_ _space_group_symop_operation_xyz\n"
                  "'x,y,z' '-x,-y,-z'\n"
                  "loop_ _atom_site_label _atom_site_type_symbol\n"
                  "_atom_site_fract_x _atom_site_fract_y _atom_site_fract_z\n"
                  "O1 O  0.3 0.0 0.0\n"
                  "O2 O -0.3 0.0 0.0\n");
  auto const structure = one_structure(text);
  CHECK(structure.cell.size() == 2);
  CHECK(structure.labels == std::vector<std::string>{"O1", "O1"});
}

TEST_CASE("a null type symbol falls back to the label", "[cif]") {
  auto const text =
      cubic_block(
                  "loop_ _atom_site_label _atom_site_type_symbol\n"
                  "_atom_site_fract_x _atom_site_fract_y _atom_site_fract_z\n"
                  "Cl1 ? 0.0 0.0 0.0\n"
                  "D1  . 0.5 0.0 0.0\n");
  auto const structure = one_structure(text);
  CHECK(structure.cell.types() ==
        Types{*data::atomic_number("Cl"), *data::atomic_number("H")});
}

TEST_CASE("the species comes from the label when no type symbol is given",
          "[cif]") {
  auto const text =
      cubic_block(
                  "loop_ _atom_site_label\n"
                  "_atom_site_fract_x _atom_site_fract_y _atom_site_fract_z\n"
                  "O1   0.0 0.0 0.0\n"
                  "Fe3+ 0.5 0.0 0.0\n"
                  "Cl1  0.0 0.5 0.0\n");
  auto const structure = one_structure(text);
  CHECK(structure.cell.types() ==
        Types{*data::atomic_number("O"), *data::atomic_number("Fe"),
              *data::atomic_number("Cl")});
}

TEST_CASE("esd-bearing numbers are read as their value", "[cif]") {
  auto const text =
      std::string{"data_x\n"
                  "_cell_length_a 4.0000(12)\n_cell_length_b 4.0\n"
                  "_cell_length_c 4.0\n"
                  "_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                  "_cell_angle_gamma 90\n"
                  "loop_ _atom_site_label\n"
                  "_atom_site_fract_x _atom_site_fract_y _atom_site_fract_z\n"
                  "Na1 0.2500( 0.0 0.0\n"};
  auto const structure = one_structure(text);
  CHECK(structure.cell.lattice().matrix()(0, 0) == Catch::Approx(4.0));
  CHECK(structure.cell.position(0)[0] == Catch::Approx(0.25));
}

TEST_CASE("the reader reports what it cannot make sense of", "[cif]") {
  auto const with = [](std::string_view body) {
    return std::string{"data_x\n"} + std::string{body};
  };
  // A cell parameter the reader needs and the block does not carry.
  CHECK(errored([&] {
    return io::read_cif(with("_cell_length_a 4.0\n_cell_length_c 4.0\n"
                             "_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                             "_cell_angle_gamma 90\n"
                             "loop_ _atom_site_label _atom_site_fract_x "
                             "_atom_site_fract_y _atom_site_fract_z\n"
                             "Na1 0 0 0\n"));
  }));
  // An element symbol no table carries.
  CHECK(errored([&] {
    return io::read_cif(with("_cell_length_a 4.0\n_cell_length_b 4.0\n"
                             "_cell_length_c 4.0\n"
                             "_cell_angle_alpha 90\n_cell_angle_beta 90\n"
                             "_cell_angle_gamma 90\n"
                             "loop_ _atom_site_label _atom_site_fract_x "
                             "_atom_site_fract_y _atom_site_fract_z\n"
                             "Zz1 0 0 0\n"));
  }));
  // A space-group symbol no setting carries.
  CHECK(errored([&] {
    return io::read_cif(nacl_block("_space_group_name_H-M_alt 'Q q q'\n"));
  }));
}

TEST_CASE("a block without a cell is skipped, not errored", "[cif]") {
  auto const structures =
      must(io::read_cif(std::string{"data_notes\n_comment 'no cell here'\n"} +
                        nacl_block(kFccSymops)));
  REQUIRE(structures.size() == 1);
  CHECK(structures.front().name == "nacl");
}

TEST_CASE("symbol lookup accepts the spellings files use", "[cif]") {
  using seitz::GroupFamily;
  auto const number = [](std::string_view symbol) {
    auto const hall = data::hall_from_hm_symbol<GroupFamily::space>(symbol);
    return hall ? data::spacegroup_type(*hall).number : 0;
  };
  CHECK(number("P 21/c") == 14);
  CHECK(number("P 1 21/c 1") == 14);
  CHECK(number("P2_1/c") == 14);
  CHECK(number("P2(1)/c") == 14);
  CHECK(number("Pm3m") == 221);
  CHECK(number("Pnma") == 62);
  CHECK(number("Pbnm") == 62);
  CHECK(number("Xyz") == 0);
  auto const hall = data::hall_from_hall_symbol<GroupFamily::space>("-P 2ybc");
  REQUIRE(hall.has_value());
  CHECK(data::spacegroup_type(*hall).number == 14);
}

// ---- the writer ------------------------------------------------------------

namespace {

// Real atomic numbers as types. Not test::rocksalt_motif, whose 0 and 1 are
// labels: the writer spells 0 as `X` and that does not read back.
[[nodiscard]] Cell nacl_cell() {
  Positions positions(2, 3);
  positions.row(0) << 0.0, 0.0, 0.0;
  positions.row(1) << 0.5, 0.5, 0.5;
  return Cell{Lattice{Matrix3d::Identity() * 5.64}, positions,
              Types{*data::atomic_number("Na"), *data::atomic_number("Cl")}};
}

} // namespace

TEST_CASE("a P1 document round-trips every atom", "[cif]") {
  auto const cell = nacl_cell();
  auto const text = io::write_cif(cell, "rocksalt");
  auto const back = one_structure(text);
  CHECK(back.name == "rocksalt");
  CHECK(back.cell.size() == cell.size());
  CHECK(back.cell.lattice().matrix().isApprox(cell.lattice().matrix(), 1e-6));
  for (auto const [row, atom] : std::views::enumerate(cell.atoms())) {
    auto const &[position, type] = atom;
    // enumerate over a view of an iota gives a difference_type of __int128,
    // which Catch cannot stream; the cell indexes by Index anyway.
    auto const i = static_cast<Index>(row);
    INFO(i);
    CHECK(back.cell.position(i).isApprox(position, 1e-5));
  }
}

TEST_CASE("a symmetrized document names the group it was determined as",
          "[cif]") {
  auto const cell = nacl_cell();
  auto const analyzer =
      analysis::SymmetryAnalyzer::from_cell(cell, io::kCifTolerance);
  auto const text = must(io::write_cif(analyzer, "rocksalt"));

  // Header states the setting three ways plus the loop; all must agree.
  auto const block = one_block(text);
  CHECK(block.value("_space_group_it_number") == "221");
  REQUIRE(block.column("_space_group_symop_operation_xyz").has_value());
  CHECK(block.column("_atom_site_wyckoff_symbol").has_value());

  auto const back = one_structure(text);
  REQUIRE(back.hall.has_value());
  CHECK(data::spacegroup_type(*back.hall).number == 221);
  CHECK(back.cell.size() == must(analyzer.standardized_cell()).size());
}

TEST_CASE("as_cif prints what write_cif returns", "[cif]") {
  auto const cell = nacl_cell();
  CHECK(std::format("{}", io::as_cif(cell, "rocksalt")) ==
        io::write_cif(cell, "rocksalt"));

  auto const analyzer =
      analysis::SymmetryAnalyzer::from_cell(cell, io::kCifTolerance);
  CHECK(std::format("{}", io::as_cif(analyzer, "rocksalt")) ==
        must(io::write_cif(analyzer, "rocksalt")));
}

TEST_CASE("a type that is not an element is written as X", "[cif]") {
  // A type outside 1..96 has no symbol, so `X`, and does not read back. Type 1
  // is a real atomic number and comes out as hydrogen.
  auto const text = io::write_cif(test::rocksalt_motif(5.64));
  auto const block = one_block(text);
  CHECK(rows(block, "_atom_site_type_symbol") ==
        std::vector<std::string>{"X", "H"});
  CHECK(errored([&] { return io::read_cif(text); }));
}

TEST_CASE("a failed determination formats as a comment, not a structure",
          "[cif]") {
  // An empty cell has nothing to determine.
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(
      Cell{Lattice{Matrix3d::Identity() * 4.0}, Positions{0, 3}, Types{}},
      io::kCifTolerance);
  auto const text = std::format("{}", io::as_cif(analyzer));
  CHECK(text.starts_with("#"));
  CHECK(must(io::parse_cif(text)).empty());
}
