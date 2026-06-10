#pragma once

#include <spglib/core/types.hpp>

#include <string_view>

// Decoder for the site-symmetry / Wyckoff database (generated tables in
// data/sitesym_tables.hpp). Port of sitesym_database.c's ssmdb_* accessors,
// 3D space-group path only (Hall numbers 1..530). The encoded tables and the
// compile-time decode are an implementation detail of sitesym_database.cpp;
// this header exposes only the decoded, Eigen-valued API.
namespace spglib::data {

// A Wyckoff position's representative coordinate operation x' = rot.x + trans
// (used to map an atom onto the canonical Wyckoff coordinate) and its
// multiplicity. Port of ssmdb_get_coordinate.
struct WyckoffCoordinate {
  Matrix3i rotation;
  Vector3d translation;
  int multiplicity;
};

[[nodiscard]] WyckoffCoordinate wyckoff_coordinate(int index);

// The half-open range [start, start + count) of global Wyckoff-position indices
// belonging to a Hall number. Port of ssmdb_get_wyckoff_indices (index > 0).
struct WyckoffRange {
  int start;
  int count;
};

[[nodiscard]] WyckoffRange wyckoff_indices(int hall_number);

// Site-symmetry symbol of a Wyckoff position (ssmdb_get_site_symmetry_symbol).
// The trailing space padding is removed at generation time, so the view points
// straight into static storage.
[[nodiscard]] std::string_view site_symmetry_symbol(int index);

} // namespace spglib::data
