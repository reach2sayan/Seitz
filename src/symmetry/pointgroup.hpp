#pragma once

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <array>
#include <optional>
#include <string_view>
#include <span>

namespace cppcrystal::symmetry {

// The 32 point groups: for each, the histogram of rotation types that
// identifies it, plus its symbols and crystal-system classification. Entry 0 is
// the empty "no point group" row.
struct PgEntry {
  std::array<int, 10> table;
  std::string_view symbol;
  std::string_view schoenflies;
  Holohedry holohedry;
  Laue laue;
};

constexpr std::array<PgEntry, 33> kPointgroupData = {{
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "", "", Holohedry::none, Laue::none},
    {{0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     "1",
     "C1",
     Holohedry::triclinic,
     Laue::laue_1},
    {{0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     "-1",
     "Ci",
     Holohedry::triclinic,
     Laue::laue_1},
    {{0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
     "2",
     "C2",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
     "m",
     "Cs",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
     "2/m",
     "C2h",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 0, 0, 1, 3, 0, 0, 0},
     "222",
     "D2",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 2, 0, 1, 1, 0, 0, 0},
     "mm2",
     "C2v",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 3, 1, 1, 3, 0, 0, 0},
     "mmm",
     "D2h",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 0, 0, 1, 1, 0, 2, 0},
     "4",
     "C4",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 2, 0, 0, 0, 1, 1, 0, 0, 0},
     "-4",
     "S4",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 2, 0, 1, 1, 1, 1, 0, 2, 0},
     "4/m",
     "C4h",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 0, 0, 0, 0, 1, 5, 0, 2, 0},
     "422",
     "D4",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 0, 0, 4, 0, 1, 1, 0, 2, 0},
     "4mm",
     "C4v",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 2, 0, 2, 0, 1, 3, 0, 0, 0},
     "-42m",
     "D2d",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 2, 0, 5, 1, 1, 5, 0, 2, 0},
     "4/mmm",
     "D4h",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 0, 0, 0, 0, 1, 0, 2, 0, 0},
     "3",
     "C3",
     Holohedry::trigonal,
     Laue::laue_3},
    {{0, 0, 2, 0, 1, 1, 0, 2, 0, 0},
     "-3",
     "C3i",
     Holohedry::trigonal,
     Laue::laue_3},
    {{0, 0, 0, 0, 0, 1, 3, 2, 0, 0},
     "32",
     "D3",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 0, 3, 0, 1, 0, 2, 0, 0},
     "3m",
     "C3v",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 2, 3, 1, 1, 3, 2, 0, 0},
     "-3m",
     "D3d",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 0, 0, 0, 1, 1, 2, 0, 2},
     "6",
     "C6",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{2, 0, 0, 1, 0, 1, 0, 2, 0, 0},
     "-6",
     "C3h",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{2, 0, 2, 1, 1, 1, 1, 2, 0, 2},
     "6/m",
     "C6h",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{0, 0, 0, 0, 0, 1, 7, 2, 0, 2},
     "622",
     "D6",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{0, 0, 0, 6, 0, 1, 1, 2, 0, 2},
     "6mm",
     "C6v",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{2, 0, 0, 4, 0, 1, 3, 2, 0, 0},
     "-6m2",
     "D3h",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{2, 0, 2, 7, 1, 1, 7, 2, 0, 2},
     "6/mmm",
     "D6h",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{0, 0, 0, 0, 0, 1, 3, 8, 0, 0},
     "23",
     "T",
     Holohedry::cubic,
     Laue::laue_m3},
    {{0, 0, 8, 3, 1, 1, 3, 8, 0, 0},
     "m-3",
     "Th",
     Holohedry::cubic,
     Laue::laue_m3},
    {{0, 0, 0, 0, 0, 1, 9, 8, 6, 0},
     "432",
     "O",
     Holohedry::cubic,
     Laue::laue_m3m},
    {{0, 6, 0, 6, 0, 1, 3, 8, 0, 0},
     "-43m",
     "Td",
     Holohedry::cubic,
     Laue::laue_m3m},
    {{0, 6, 8, 9, 1, 1, 9, 8, 6, 0},
     "m-3m",
     "Oh",
     Holohedry::cubic,
     Laue::laue_m3m},
}};

// Point-group data (symbol, Schoenflies, holohedry, Laue) for number 1..32;
// number 0 / out of range gives an empty PointGroup.
[[nodiscard]] constexpr PointGroup pointgroup_by_number(int number) noexcept {
  if (number < 1 || number > 32) {
    return {};
  }
  PgEntry const &e = kPointgroupData[static_cast<std::size_t>(number)];
  // CrystalClass enumerators are aligned to the point-group numbering (1..32).
  return {number, e.symbol, e.schoenflies, e.holohedry, e.laue,
          static_cast<CrystalClass>(number)};
}

// The 32 crystal classes, in the International Tables order. Pinned at compile
// time so a drift in the table breaks the build rather than a test.
static_assert([] {
  constexpr std::array<std::string_view, 33> kSymbols = {
      "",      "1",   "-1",  "2",     "m",    "2/m",  "222",   "mm2",
      "mmm",   "4",   "-4",  "4/m",   "422",  "4mm",  "-42m",  "4/mmm",
      "3",     "-3",  "32",  "3m",    "-3m",  "6",    "-6",    "6/m",
      "622",   "6mm", "-6m2", "6/mmm", "23",  "m-3",  "432",   "-43m",
      "m-3m"};
  for (int n = 1; n <= 32; ++n) {
    if (pointgroup_by_number(n).symbol !=
        kSymbols[static_cast<std::size_t>(n)]) {
      return false;
    }
  }
  return pointgroup_by_number(0).number == 0 &&
         pointgroup_by_number(33).number == 0;
}());

// Point group together with the integer change-of-basis matrix that brings the
// rotations into the conventional setting (columns are the chosen axes). With
// `aperiodic_axis` set (a layer cell), the cubic point groups are rejected and
// the conventional axes are chosen so the aperiodic axis becomes c.
struct PointgroupTransform {
  PointGroup pointgroup;
  Matrix3i transformation{Matrix3i::Zero()};
};

// The family is a compile-time parameter: only the layer path rejects the
// cubic point groups and sorts the axes so the aperiodic one becomes c.
// `layer_axis` is still data — a layer cell's aperiodic axis is whichever of
// the three the input basis put it on — and is ignored for GroupFamily::space.
template <GroupFamily F>
[[nodiscard]] Result<PointgroupTransform>
identify_point_group(std::span<Matrix3i const> rotations,
                     std::optional<int> layer_axis = std::nullopt);

extern template Result<PointgroupTransform>
identify_point_group<GroupFamily::space>(std::span<Matrix3i const>,
                                         std::optional<int>);
extern template Result<PointgroupTransform>
identify_point_group<GroupFamily::layer>(std::span<Matrix3i const>,
                                         std::optional<int>);

} // namespace cppcrystal::symmetry
