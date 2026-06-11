#pragma once

#include <spglib/core/point_group.hpp> // Holohedry

#include <array>
#include <cstddef>

// Compile-time classification of Hall numbers (1..530). Every property here is a
// pure integer function of the Hall number over a fixed, finite domain, so it is
// resolved entirely at compile time into a single lookup table (kHallClass),
// following the house constexpr-IIFE idiom used by kCatalog (spg_database.hpp)
// and kDecodedOps (spg_database.cpp). This is the one source of truth for the
// crystal-system bucket and the rhombohedral subsets that hall_symbol.cpp and
// spacegroup.cpp previously hardcoded in several places. No Eigen here — the
// Eigen-valued symmetry operations stay in the lazily-built runtime cache.
namespace spglib::data {

// Pure integer-derived properties of a Hall number (1..530).
struct HallClass {
  Holohedry system = Holohedry::none;
  bool rhombohedral = false;       // R-centered subset of the trigonal range
  bool rhombo_hex_setting = false; // hP (hexagonal) setting of the rhombo subset
};

// Crystal system from the Hall-number range (hall_symbol.c dispatch ranges). The
// seven non-`none` Holohedry values map exactly onto these ranges.
[[nodiscard]] constexpr Holohedry hall_crystal_system(int h) noexcept {
  if (489 <= h && h <= 530)
    return Holohedry::cubic;
  if (462 <= h && h <= 488)
    return Holohedry::hexagonal;
  if (430 <= h && h <= 461)
    return Holohedry::trigonal;
  if (349 <= h && h <= 429)
    return Holohedry::tetragonal;
  if (108 <= h && h <= 348)
    return Holohedry::orthorhombic;
  if (3 <= h && h <= 107)
    return Holohedry::monoclinic;
  if (1 <= h && h <= 2)
    return Holohedry::triclinic;
  return Holohedry::none;
}

// Crystal system from a point-group number (1..32), by the standard contiguous
// ranges. Used to classify layer-group Hall settings (negative hall numbers),
// which carry no 3D hall-range bucket but do carry a point-group number.
[[nodiscard]] constexpr Holohedry holohedry_from_pointgroup(int pg) noexcept {
  if (pg < 1 || pg > 32)
    return Holohedry::none;
  if (pg <= 2)
    return Holohedry::triclinic;
  if (pg <= 5)
    return Holohedry::monoclinic;
  if (pg <= 8)
    return Holohedry::orthorhombic;
  if (pg <= 15)
    return Holohedry::tetragonal;
  if (pg <= 20)
    return Holohedry::trigonal;
  if (pg <= 27)
    return Holohedry::hexagonal;
  return Holohedry::cubic;
}

// R-centered (rhombohedral) subset of the trigonal range (hall_symbol.c).
[[nodiscard]] constexpr bool is_rhombohedral_hall(int h) noexcept {
  switch (h) {
  case 433:
  case 434:
  case 436:
  case 437:
  case 444:
  case 445:
  case 450:
  case 451:
  case 452:
  case 453:
  case 458:
  case 459:
  case 460:
  case 461:
    return true;
  default:
    return false;
  }
}

// The hexagonal (hP) setting within the rhombohedral subset (hall_symbol.c's
// is_rhombo_h_setting / spacegroup.c's match_db_rhombo `hex` test).
[[nodiscard]] constexpr bool is_rhombo_hex_setting(int h) noexcept {
  switch (h) {
  case 433:
  case 436:
  case 444:
  case 450:
  case 452:
  case 458:
  case 460:
    return true;
  default:
    return false;
  }
}

// All 530 Hall numbers classified once at compile time. Index 0 is the
// out-of-range sentinel (Holohedry::none, both flags false).
inline constexpr std::array<HallClass, 531> kHallClass = [] {
  std::array<HallClass, 531> t{};
  for (int h = 1; h <= 530; ++h)
    t[static_cast<std::size_t>(h)] =
        HallClass{hall_crystal_system(h), is_rhombohedral_hall(h),
                  is_rhombo_hex_setting(h)};
  return t;
}();

[[nodiscard]] constexpr HallClass const &hall_class(int h) noexcept {
  bool const in_range = h >= 1 && h <= 530;
  return kHallClass[static_cast<std::size_t>(in_range ? h : 0)];
}

// Compile-time guards locking the consolidated sets to their reference shape;
// these break the build if the three properties ever drift out of agreement.
namespace detail {
constexpr int count_rhombohedral = [] {
  int n = 0;
  for (auto const &c : kHallClass)
    n += c.rhombohedral ? 1 : 0;
  return n;
}();
constexpr int count_rhombo_hex = [] {
  int n = 0;
  for (auto const &c : kHallClass)
    n += c.rhombo_hex_setting ? 1 : 0;
  return n;
}();
constexpr bool flags_imply_trigonal = [] {
  for (auto const &c : kHallClass) {
    // rhombo_hex ⟹ rhombohedral ⟹ trigonal
    if (c.rhombo_hex_setting && !c.rhombohedral)
      return false;
    if (c.rhombohedral && c.system != Holohedry::trigonal)
      return false;
  }
  return true;
}();
} // namespace detail

static_assert(detail::count_rhombohedral == 14);
static_assert(detail::count_rhombo_hex == 7);
static_assert(detail::flags_imply_trigonal);
static_assert(hall_class(0).system == Holohedry::none);   // sentinel
static_assert(hall_class(531).system == Holohedry::none); // out of range
static_assert(hall_class(1).system == Holohedry::triclinic);
static_assert(hall_class(530).system == Holohedry::cubic);
static_assert(hall_class(433).rhombohedral && hall_class(433).rhombo_hex_setting);
static_assert(hall_class(434).rhombohedral && !hall_class(434).rhombo_hex_setting);

} // namespace spglib::data
