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

// A scalar item is a one-row column. std::less<> for heterogeneous lookup.
using CifColumns =
    boost::container::flat_map<std::string, std::vector<std::string>,
                               std::less<>>;

// One `data_` block as written: tags lowercased, values uninterpreted.
struct CifBlock {
  std::string name;
  CifColumns columns;
  std::vector<std::vector<std::string>> loops; // which tags shared a loop_

  // nullopt when the block does not carry `tag`.
  [[nodiscard]] std::optional<std::span<std::string const>>
  column(std::string_view tag) const;

  // First row of `tag`; nullopt for absent and for CIF's nulls '?' and '.'.
  [[nodiscard]] std::optional<std::string_view> value(std::string_view tag) const;
};

// The blocks of a CIF 1.1 document, in file order. Errors e_cif_syntax with
// the 1-based line and column the parse stopped at.
[[nodiscard]] Result<std::vector<CifBlock>> parse_cif(std::string_view text);

// A partially occupied or mixed site, collapsed to its majority species.
// Reported rather than silently applied: the structure is not the one on disk.
struct OccupancyCollapse {
  std::string kept;                  // the label that survived
  double occupancy = 1.0;            // its stated occupancy
  std::vector<std::string> dropped;  // the labels that shared the site
};

// One block read as crystallography.
struct CifStructure {
  std::string name;
  Cell cell;
  // nullopt when the file named no setting (P1 assumed), and when a symop
  // loop matched none -- the operations still expanded the cell.
  std::optional<HallNumber> hall;
  std::vector<std::string> labels; // per cell atom, its asymmetric-unit label
  std::vector<OccupancyCollapse> collapsed;
};

// Five-decimal coordinates put two images of a site ~1e-4 apart, which the
// 1e-5 search default would call different atoms. Do not substitute Tolerance{}.
inline constexpr Tolerance kCifTolerance{.symprec = 1e-3,
                                        .angle_tolerance = std::nullopt};

// One block as a structure, symmetry expanded over the asymmetric unit. Errors
// e_cif_missing, e_unknown_element, e_unknown_spacegroup_symbol,
// e_invalid_xyz, e_invalid_lattice, e_empty_cell.
[[nodiscard]] Result<CifStructure> structure_of(CifBlock const &block,
                                                Tolerance tol = kCifTolerance);

// Every block carrying a cell. Blocks without `_cell_length_a` are skipped,
// not errored: real files carry commentary blocks beside their structures.
[[nodiscard]] Result<std::vector<CifStructure>>
read_cif(std::string_view text, Tolerance tol = kCifTolerance);

// The cell in P1, every atom written out. Total. A type that is not a
// tabulated element is written `X` and does not survive a round trip.
[[nodiscard]] std::string write_cif(Cell const &cell,
                                    std::string_view name = "seitz");

// The standardized cell, its setting's database operations, and one atom per
// orbit with Wyckoff letter and multiplicity. Errors as the determination does.
[[nodiscard]] Result<std::string>
write_cif(analysis::SymmetryAnalyzer const &analyzer,
          std::string_view name = "seitz");

// The same rendering through std::format/std::print: `std::print("{}",
// as_cif(cell))`. A view -- `subject` must outlive the format call.
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

// Formatting cannot fail and the determination can, so a failure writes one
// CIF comment -- valid CIF parsing back to zero blocks, never a structure.
// For the error tag itself, call write_cif and handle the Result.
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
