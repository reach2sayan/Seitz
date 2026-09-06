#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/core/cell.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/io/cif.hpp>

#include "errors.hpp" // unwrap

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Text in, text out: no std::filesystem crosses this boundary -- the
// pure-Python wrapper reads a Path far more naturally than a binding could.
namespace seitz::python {

void bind_io(py::module_ &m) {
  py::class_<io::CifBlock>(
      m, "CifBlock",
      "One `data_` block as written: tags lowercased, values still the file's "
      "own strings. Nothing here is interpreted as crystallography.")
      .def_readonly("name", &io::CifBlock::name)
      .def_property_readonly(
          "columns",
          [](io::CifBlock const &self) {
            // flat_map has no caster, and a dict is what a caller wants.
            py::dict out;
            for (auto const &[tag, rows] : self.columns) {
              out[py::str(tag)] = py::cast(rows);
            }
            return out;
          },
          "Every tag's rows, as a dict. A scalar item is a one-row column.")
      .def_readonly("loops", &io::CifBlock::loops,
                    "The tags of each `loop_`, in file order.")
      .def(
          "column",
          [](io::CifBlock const &self, std::string_view tag)
              -> std::optional<std::vector<std::string>> {
            auto const rows = self.column(tag);
            if (!rows) {
              return std::nullopt;
            }
            return std::vector<std::string>{rows->begin(), rows->end()};
          },
          py::arg("tag"), "The rows of `tag`, or None if the block lacks it.")
      .def(
          "value",
          [](io::CifBlock const &self,
             std::string_view tag) -> std::optional<std::string> {
            auto const text = self.value(tag);
            if (!text) {
              return std::nullopt;
            }
            return std::string{*text};
          },
          py::arg("tag"),
          "The scalar reading of `tag`: its first row, with CIF's '?' and '.' "
          "reported as None like an absent tag.")
      .def("__repr__", [](io::CifBlock const &self) {
        return "CifBlock(name=" + self.name +
               ", tags=" + std::to_string(self.columns.size()) + ")";
      });

  py::class_<io::OccupancyCollapse>(
      m, "OccupancyCollapse",
      "A site the file could not state exactly: a partially occupied or "
      "mixed-species position, collapsed to its majority species.")
      .def_readonly("kept", &io::OccupancyCollapse::kept,
                    "The label that survived.")
      .def_readonly("occupancy", &io::OccupancyCollapse::occupancy,
                    "Its stated occupancy.")
      .def_readonly("dropped", &io::OccupancyCollapse::dropped,
                    "The labels that shared the site.")
      .def("__repr__", [](io::OccupancyCollapse const &self) {
        return "OccupancyCollapse(kept=" + self.kept + ", dropped=" +
               std::to_string(self.dropped.size()) + ")";
      });

  py::class_<io::CifStructure>(m, "CifStructure",
                               "One block read as crystallography.")
      .def_readonly("name", &io::CifStructure::name)
      .def_readonly("cell", &io::CifStructure::cell)
      .def_readonly("hall", &io::CifStructure::hall,
                    "The setting the file named, or None when it named none "
                    "and P1 was assumed.")
      .def_readonly("labels", &io::CifStructure::labels,
                    "Per cell atom, its asymmetric-unit label.")
      .def_readonly("collapsed", &io::CifStructure::collapsed,
                    "What the reader had to decide; empty for an ordered file.")
      .def("__repr__", [](io::CifStructure const &self) {
        return "CifStructure(name=" + self.name +
               ", atoms=" + std::to_string(self.cell.size()) + ")";
      });

  m.def(
      "parse_cif",
      [](std::string_view text) {
        return unwrap([&] { return io::parse_cif(text); });
      },
      py::arg("text"),
      "The blocks of a CIF document, in file order, uninterpreted.");

  m.def(
      "read_cif",
      [](std::string_view text, Tolerance const &tol) {
        return unwrap([&] { return io::read_cif(text, tol); });
      },
      py::arg("text"), py::arg("tolerance") = io::kCifTolerance,
      "Every block that carries a cell, read as a structure.");

  m.def(
      "write_cif_cell",
      [](Cell const &cell, std::string_view name) {
        return io::write_cif(cell, name);
      },
      py::arg("cell"), py::arg("name"),
      "The cell as a CIF document, in P1.");

  m.def(
      "write_cif_analyzer",
      [](analysis::SymmetryAnalyzer const &analyzer, std::string_view name) {
        return unwrap([&] { return io::write_cif(analyzer, name); });
      },
      py::arg("analyzer"), py::arg("name"),
      "The determination as a CIF document: standardized cell, database "
      "operations, one representative per orbit.");

  m.attr("CIF_SYMPREC") = io::kCifTolerance.symprec;
}

} // namespace seitz::python
