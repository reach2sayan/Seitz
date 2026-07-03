#pragma once

#include <cppcrystal/core/types.hpp>

#include <boost/container/static_vector.hpp>

#include <string_view>

namespace cppcrystal {

// Crystal system / holohedry.
enum class Holohedry {
  none,
  triclinic,
  monoclinic,
  orthorhombic,
  tetragonal,
  trigonal,
  hexagonal,
  cubic,
};

// Laue class.
enum class Laue {
  none,
  laue_1,    // -1
  laue_2m,   // 2/m
  laue_mmm,  // mmm
  laue_4m,   // 4/m
  laue_4mmm, // 4/mmm
  laue_3,    // -3
  laue_3m,   // -3m
  laue_6m,   // 6/m
  laue_6mmm, // 6/mmm
  laue_m3,   // m-3
  laue_m3m,  // m-3m
};

// The 32 crystallographic point groups (crystal classes), named by their
// Schoenflies symbol. Enumerator values match the international point-group
// numbering (1..32, 0 = none) used by pointgroup_by_number /
// identify_pointgroup_number, so static_cast<CrystalClass>(number) round-trips.
enum class CrystalClass {
  none,
  c1,  // 1
  ci,  // -1
  c2,  // 2
  cs,  // m
  c2h, // 2/m
  d2,  // 222
  c2v, // mm2
  d2h, // mmm
  c4,  // 4
  s4,  // -4
  c4h, // 4/m
  d4,  // 422
  c4v, // 4mm
  d2d, // -42m
  d4h, // 4/mmm
  c3,  // 3
  c3i, // -3
  d3,  // 32
  c3v, // 3m
  d3d, // -3m
  c6,  // 6
  c3h, // -6
  c6h, // 6/m
  d6,  // 622
  c6v, // 6mm
  d3h, // -6m2
  d6h, // 6/mmm
  t,   // 23
  th,  // m-3
  o,   // 432
  td,  // -43m
  oh,  // m-3m
};

struct PointGroup {
  int number = 0;               // 1..32; 0 means unknown / not found
  std::string_view symbol;      // international (Hermann-Mauguin) symbol
  std::string_view schoenflies; // Schoenflies symbol
  Holohedry holohedry = Holohedry::none;
  Laue laue = Laue::none;
  CrystalClass crystal_class = CrystalClass::none;
};

using PointSymmetry = boost::container::static_vector<Matrix3i, 48>;

} // namespace cppcrystal
