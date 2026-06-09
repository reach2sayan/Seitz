#pragma once

#include <spglib/core/types.hpp>

#include <boost/container/static_vector.hpp>

#include <string_view>

namespace spglib {

// Crystal system / holohedry (spglib Holohedry enum).
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

// Laue class (spglib Laue enum).
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

struct PointGroup {
  int number = 0;               // 1..32; 0 means unknown / not found
  std::string_view symbol;      // international (Hermann-Mauguin) symbol
  std::string_view schoenflies; // Schoenflies symbol
  Holohedry holohedry = Holohedry::none;
  Laue laue = Laue::none;
};

using PointSymmetry = boost::container::static_vector<Matrix3i, 48>;

} // namespace spglib
