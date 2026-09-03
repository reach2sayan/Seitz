#pragma once

// Loader for spglib's reference test corpus: one input cell per space group at
// <ref>/test/functional/python/data/<system>/unitcell_<N>.yaml. The path comes
// from the SPGLIB_REF_DATA_DIR compile definition (set in tests/CMakeLists.txt
// from the FetchContent source dir). A tiny hand-rolled parser is enough — the
// files are machine-generated and regular (lattice = 3 bracketed float triples,
// points = number + coordinates, plus the expected space_group.number).

#include <cppcrystal/core/cell.hpp>

#include <cctype>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <system_error>
#include <vector>

namespace cppcrystal::oracle {

struct CorpusEntry {
  Cell cell;
  int space_group_number = 0;
  std::string name;
};

// Every numeric token on a line (handles signs, decimals, exponents). The YAML
// list dash "- [" yields a lone "-" which fails to parse and is skipped.
inline std::vector<double> extract_numbers(std::string const &s) {
  std::vector<double> out;
  std::size_t i = 0;
  while (i < s.size()) {
    char const c = s[i];
    bool const starts =
        c == '-' || c == '+' || c == '.' ||
        std::isdigit(static_cast<unsigned char>(c)) != 0;
    if (!starts) {
      ++i;
      continue;
    }
    std::size_t j = i;
    while (j < s.size()) {
      char const d = s[j];
      if (std::isdigit(static_cast<unsigned char>(d)) != 0 || d == '.' ||
          d == '-' || d == '+' || d == 'e' || d == 'E') {
        ++j;
      } else {
        break;
      }
    }
    try {
      out.push_back(std::stod(s.substr(i, j - i)));
    } catch (...) { // lone "-", "+", "." etc.
    }
    i = j;
  }
  return out;
}

inline std::optional<CorpusEntry>
parse_unitcell(std::filesystem::path const &path) {
  std::ifstream file(path);
  if (!file) {
    return std::nullopt;
  }

  int space_group = 0;
  std::vector<Vector3d> lattice_rows;
  std::vector<int> types;
  std::vector<Vector3d> coordinates;

  enum class Section { none, space_group, lattice, points } section = Section::none;
  for (std::string line; std::getline(file, line);) {
    auto const has = [&](char const *key) {
      return line.find(key) != std::string::npos;
    };
    if (has("space_group:")) {
      section = Section::space_group;
      continue;
    }
    if (has("lattice:")) {
      section = Section::lattice;
      continue;
    }
    if (has("points:")) {
      section = Section::points;
      continue;
    }

    switch (section) {
    case Section::space_group:
      if (has("number:")) {
        auto const n = extract_numbers(line);
        if (!n.empty()) {
          space_group = static_cast<int>(n[0]);
        }
        section = Section::none;
      }
      break;
    case Section::lattice:
      if (has("[")) {
        auto const n = extract_numbers(line);
        if (n.size() >= 3) {
          lattice_rows.emplace_back(n[0], n[1], n[2]);
        }
      }
      break;
    case Section::points:
      if (has("number:")) {
        auto const n = extract_numbers(line);
        if (!n.empty()) {
          types.push_back(static_cast<int>(n[0]));
        }
      } else if (has("coordinates:") && has("[")) {
        auto const n = extract_numbers(line);
        if (n.size() >= 3) {
          coordinates.emplace_back(n[0], n[1], n[2]);
        }
      }
      break;
    case Section::none:
      break;
    }
  }

  if (lattice_rows.size() != 3 || coordinates.empty() ||
      coordinates.size() != types.size()) {
    return std::nullopt;
  }

  Matrix3d lattice;
  for (int c = 0; c < 3; ++c) {
    // YAML lattice row i is basis vector i; our Matrix3d stores basis vectors as
    // columns.
    lattice.col(c) = lattice_rows[static_cast<std::size_t>(c)];
  }
  Positions positions(static_cast<Index>(coordinates.size()), 3);
  for (std::size_t a = 0; a < coordinates.size(); ++a) {
    positions.row(static_cast<Index>(a)) = coordinates[a].transpose();
  }
  return CorpusEntry{Cell(Lattice{lattice}, positions, types), space_group,
                     path.filename().string()};
}

// All reference input cells (one per space group, ~230 total).
inline std::vector<CorpusEntry> load_corpus() {
  std::vector<CorpusEntry> out;
  std::filesystem::path const base(SPGLIB_REF_DATA_DIR);
  for (char const *system :
       {"triclinic", "monoclinic", "orthorhombic", "tetragonal", "trigonal",
        "hexagonal", "cubic"}) {
    std::filesystem::path const dir = base / system;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) {
      continue;
    }
    for (auto const &entry : std::filesystem::directory_iterator(dir, ec)) {
      auto const filename = entry.path().filename().string();
      if (filename.rfind("unitcell_", 0) == 0 &&
          entry.path().extension() == ".yaml") {
        if (auto cell = parse_unitcell(entry.path())) {
          out.push_back(std::move(*cell));
        }
      }
    }
  }
  return out;
}

} // namespace cppcrystal::oracle
