#pragma once

#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/tolerance.hpp>

#include <boost/container/flat_map.hpp>

#include <format>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::io {

// A tag's values. A scalar item is a one-row column, so nothing downstream has
// to ask which of the two forms a file happened to use. std::less<> so a
// std::string_view tag looks up without allocating.
using CifColumns =
    boost::container::flat_map<std::string, std::vector<std::string>,
                               std::less<>>;

// One `data_` block of a CIF document, as written: tags lowercased (CIF tags
// are case-insensitive), values still the file's own strings. Nothing here is
// interpreted as crystallography -- that is structure_of()'s job.
struct CifBlock {
  std::string name;
  CifColumns columns;
  // The tags of each `loop_`, in file order: which columns belong together,
  // which the document model alone no longer says.
  std::vector<std::vector<std::string>> loops;

  // The rows of `tag`, or nullopt when the block does not carry it.
  [[nodiscard]] std::optional<std::span<std::string const>>
  column(std::string_view tag) const;

  // The scalar reading of `tag`: its first row, with CIF's two null tokens
  // ('?' unknown, '.' inapplicable) reported as nullopt like an absent tag.
  [[nodiscard]] std::optional<std::string_view> value(std::string_view tag) const;
};

// The blocks of a CIF 1.1 document, in file order. Permissive by design: no
// expectation points, so a file that stops making sense reports where rather
// than throwing. Errors e_cif_syntax carrying the 1-based line and column the
// parse stopped at.
[[nodiscard]] Result<std::vector<CifBlock>> parse_cif(std::string_view text);

// What the reader had to decide for a site the file could not state exactly:
// a partially occupied or mixed-species position collapses to its majority
// species. Reported rather than silently applied -- a caller refining
// occupancies needs to know the structure it got is not the one on disk.
struct OccupancyCollapse {
  std::string kept;                  // the label that survived
  double occupancy = 1.0;            // its stated occupancy
  std::vector<std::string> dropped;  // the labels that shared the site
};

// One block read as crystallography.
struct CifStructure {
  std::string name;
  Cell cell;
  // The setting the file named, or nullopt when it named none and P1 was
  // assumed. A symop loop that matches no Hall setting also leaves this unset:
  // the operations still expanded the cell.
  std::optional<HallNumber> hall;
  std::vector<std::string> labels; // per cell atom, its asymmetric-unit label
  std::vector<OccupancyCollapse> collapsed;
};

// CIF coordinates are written to five decimals, which puts two images of one
// site as much as ~1e-4 apart; the search default of 1e-5 would call them
// different atoms. This is the tolerance the reader means, and tests should
// not replace it with Tolerance{}.
inline constexpr Tolerance kCifTolerance{.symprec = 1e-3,
                                        .angle_tolerance = std::nullopt};

// One block as a structure: cell parameters, sites, and the symmetry expanded
// over the asymmetric unit. Errors e_cif_missing (a tag the reader needs),
// e_unknown_element, e_unknown_spacegroup_symbol, e_invalid_xyz,
// e_invalid_lattice, e_empty_cell.
[[nodiscard]] Result<CifStructure> structure_of(CifBlock const &block,
                                                Tolerance tol = kCifTolerance);

// Every block of a document that carries a cell, read as a structure. Blocks
// without `_cell_length_a` are skipped rather than errored: real files carry
// commentary blocks beside their structures.
[[nodiscard]] Result<std::vector<CifStructure>>
read_cif(std::string_view text, Tolerance tol = kCifTolerance);

// The cell as a CIF document, in P1: every atom written out, one symmetry
// operation, no orbit folded. Total -- a Cell always has a lattice and atoms.
// Types that are not tabulated elements are written as `X` and do not survive
// a round trip.
[[nodiscard]] std::string write_cif(Cell const &cell,
                                    std::string_view name = "seitz");

// The determination as a CIF document: the standardized conventional cell, the
// database operations of its Hall setting, and one representative atom per
// orbit with its Wyckoff letter and multiplicity. Errors whatever the
// determination errors.
[[nodiscard]] Result<std::string>
write_cif(analysis::SymmetryAnalyzer const &analyzer,
          std::string_view name = "seitz");

// The same rendering reached through std::format / std::print, so writing a
// CIF and printing one are one API:
//
//   std::string const text = io::write_cif(cell, "nacl");
//   std::print("{}", io::as_cif(cell, "nacl"));
//
// A view, not a copy: `subject` must outlive the format call, which is always
// true for the intended `std::print("{}", as_cif(x))` spelling.
template <class T> struct AsCif {
  T const &subject;
  std::string_view name;
};

[[nodiscard]] inline AsCif<Cell> as_cif(Cell const &cell,
                                        std::string_view name = "seitz") {
  return {cell, name};
}
[[nodiscard]] inline AsCif<analysis::SymmetryAnalyzer>
as_cif(analysis::SymmetryAnalyzer const &analyzer,
       std::string_view name = "seitz") {
  return {analyzer, name};
}

} // namespace seitz::io

template <>
struct std::formatter<seitz::io::AsCif<seitz::Cell>>
    : std::formatter<std::string_view> {
  auto format(seitz::io::AsCif<seitz::Cell> const &document,
              std::format_context &ctx) const {
    return std::formatter<std::string_view>::format(
        seitz::io::write_cif(document.subject, document.name), ctx);
  }
};

// Formatting cannot fail, and the determination can. A failed one writes a
// single CIF comment naming the failure: valid CIF that parses back to zero
// blocks, so it can never be mistaken for a structure. A caller that needs the
// error tag itself calls write_cif and handles the Result.
template <>
struct std::formatter<seitz::io::AsCif<seitz::analysis::SymmetryAnalyzer>>
    : std::formatter<std::string_view> {
  auto format(seitz::io::AsCif<seitz::analysis::SymmetryAnalyzer> const &document,
              std::format_context &ctx) const {
    auto text = seitz::io::write_cif(document.subject, document.name);
    return std::formatter<std::string_view>::format(
        text ? *text : std::string{"# seitz: the determination failed\n"}, ctx);
  }
};

#pragma GCC visibility pop
