#pragma once

#include "core/testable.hpp"
#include <seitz/core/keys.hpp>
#include <seitz/core/types.hpp>

#include <string_view>
#include <vector>

// Decoder for the site-symmetry / Wyckoff database (generated tables in
// src/data/sitesym_tables.hpp). The
// encoded tables and the compile-time decode are an implementation detail of
// sitesym_database.cpp; this header exposes only the decoded, Eigen-valued API.
namespace seitz::data {

// A Wyckoff position's representative coordinate operation x' = rot.x + trans
// (used to map an atom onto the canonical Wyckoff coordinate) and its
// multiplicity.
struct WyckoffCoordinate {
  Matrix3i rotation;
  Vector3d translation;
  int multiplicity;
};

[[nodiscard]] SEITZ_TESTABLE WyckoffCoordinate
wyckoff_coordinate(int index);

// The half-open range [start, start + count) of global Wyckoff-position indices
// belonging to a Hall number.
struct WyckoffRange {
  int start;
  int count;
};

[[nodiscard]] SEITZ_TESTABLE WyckoffRange wyckoff_indices(HallNumber hall);

// Site-symmetry symbol of a Wyckoff position. The trailing space padding is
// removed at generation time, so the view points straight into static storage.
[[nodiscard]] SEITZ_TESTABLE std::string_view
site_symmetry_symbol(int index);

// One Wyckoff position of a Hall setting: its global database index (the key
// for wyckoff_coordinate / site_symmetry_symbol) paired with its Wyckoff
// letter. The database lists positions from the general position downward, so
// the letter is the reverse offset within the Hall range (0 = 'a' = the most
// special position) — see wyckoff_entries.
struct WyckoffEntry {
  int global_index;
  int letter; // 0 = 'a'
};

// The Wyckoff positions of a Hall setting, ordered by ascending letter
// (a, b, c, ...). Bundles index + letter so callers iterate one list of objects
// rather than re-deriving the letter offset from a WyckoffRange. The letter
// convention (count - offset - 1) matches the one used when assigning
// Dataset::wyckoffs in refine/site_symmetry.cpp.
[[nodiscard]] std::vector<WyckoffEntry> wyckoff_entries(HallNumber hall);

} // namespace seitz::data
