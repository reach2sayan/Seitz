#pragma once

#include <cppcrystal/core/point_group.hpp> // Holohedry

#include <array>
#include <cstddef>

namespace cppcrystal::data {

// Pure integer-derived properties of a Hall number (1..530).
struct HallClass {
  Holohedry system = Holohedry::none;
  bool rhombohedral = false; // R-centered subset of the trigonal range
  bool rhombo_hex_setting =
      false; // hP (hexagonal) setting of the rhombo subset
};

[[nodiscard]] constexpr Holohedry hall_crystal_system(int h) noexcept {
  if (489 <= h && h <= 530) {
    return Holohedry::cubic;
  } else if (462 <= h && h <= 488) {
    return Holohedry::hexagonal;
  } else if (430 <= h && h <= 461) {
    return Holohedry::trigonal;
  } else if (349 <= h && h <= 429) {
    return Holohedry::tetragonal;
  } else if (108 <= h && h <= 348) {
    return Holohedry::orthorhombic;
  } else if (3 <= h && h <= 107) {
    return Holohedry::monoclinic;
  } else if (1 <= h && h <= 2) {
    return Holohedry::triclinic;
  }
  return Holohedry::none;
}

[[nodiscard]] constexpr Holohedry holohedry_from_pointgroup(int pg) noexcept {
  if (pg < 1 || pg > 32) {
    return Holohedry::none;
  } else if (pg <= 2) {
    return Holohedry::triclinic;
  } else if (pg <= 5) {
    return Holohedry::monoclinic;
  } else if (pg <= 8) {
    return Holohedry::orthorhombic;
  } else if (pg <= 15) {
    return Holohedry::tetragonal;
  } else if (pg <= 20) {
    return Holohedry::trigonal;
  } else if (pg <= 27) {
    return Holohedry::hexagonal;
  }
  return Holohedry::cubic;
}

// R-centered (rhombohedral) subset of the trigonal range
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

// The hexagonal (hP) setting within the rhombohedral subset
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

// All 530 Hall numbers classified once at compile time.
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

namespace detail {
constexpr int count_rhombohedral = [] {
  return std::ranges::count(kHallClass, true, &HallClass::rhombohedral);
}();
constexpr int count_rhombo_hex = [] {
  return std::ranges::count_if(kHallClass, &HallClass::rhombo_hex_setting);
}();
constexpr bool flags_imply_trigonal = [] {
  return std::ranges::all_of(kHallClass, [](auto const &c) {
    // rhombo_hex ⟹ rhombohedral ⟹ trigonal
    return (!c.rhombo_hex_setting || c.rhombohedral) &&
           (!c.rhombohedral || c.system == Holohedry::trigonal);
  });
}();
} // namespace detail

// clang-format off
static_assert(detail::count_rhombohedral == 14);
static_assert(detail::count_rhombo_hex == 7);
static_assert(detail::flags_imply_trigonal);
static_assert(hall_class(0).system == Holohedry::none);   // sentinel
static_assert(hall_class(531).system == Holohedry::none); // out of range
static_assert(hall_class(1).system == Holohedry::triclinic);
static_assert(hall_class(530).system == Holohedry::cubic);
static_assert(hall_class(433).rhombohedral && hall_class(433).rhombo_hex_setting);
static_assert(hall_class(434).rhombohedral && !hall_class(434).rhombo_hex_setting);
// clang-format on
} // namespace cppcrystal::data
