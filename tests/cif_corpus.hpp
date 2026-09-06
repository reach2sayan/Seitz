#pragma once

// PyXtal's two CIF corpora: <pyxtal>/miscellaneous/cifs (209 files named after
// their space group, `=` standing in for `/`) and <pyxtal>/database/cifs (77
// real structures among non-CIF files). Path from SEITZ_PYXTAL_CIF_DIR.
//
// The reader takes text, not a path, so opening the file is the suite's job.

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

// Every `.cif` in one corpus directory, name-ordered for reproducibility.
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
