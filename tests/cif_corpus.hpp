#pragma once

// Loader for PyXtal's two CIF corpora, at <pyxtal>/miscellaneous/cifs (209
// files named after the space group they were built for, `P6_3=mmc.cif` with
// `=` standing in for `/`) and <pyxtal>/database/cifs (77 real structures plus
// some non-CIF files). The path comes from the SEITZ_PYXTAL_CIF_DIR compile
// definition, set in tests/CMakeLists.txt from the FetchContent source dir.
//
// The reader takes text, not a path -- there is no std::filesystem in the
// public API -- so reading the file is the test suite's job, and it is these
// five lines.

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace seitz::oracle {

struct CifFile {
  std::string name; // stem, `/` restored from `=`
  std::string text;
};

[[nodiscard]] inline std::string read_file(std::filesystem::path const &path) {
  std::ifstream in(path, std::ios::binary);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return std::move(buffer).str();
}

// Every `.cif` in one corpus directory, in name order so a failure is
// reproducible.
[[nodiscard]] inline std::vector<CifFile> cif_corpus(std::string_view corpus) {
  std::filesystem::path const directory =
      std::filesystem::path{SEITZ_PYXTAL_CIF_DIR} / corpus / "cifs";
  std::vector<CifFile> files;
  for (auto const &entry : std::filesystem::directory_iterator{directory}) {
    if (entry.path().extension() != ".cif") {
      continue;
    }
    std::string name = entry.path().stem().string();
    std::ranges::replace(name, '=', '/');
    files.push_back(CifFile{.name = std::move(name),
                            .text = read_file(entry.path())});
  }
  std::ranges::sort(files, {}, &CifFile::name);
  return files;
}

} // namespace seitz::oracle
