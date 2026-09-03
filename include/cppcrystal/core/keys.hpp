#pragma once

namespace cppcrystal {

// The two crystallographic families the determination pipeline handles. A
// layer group is 2D-periodic with its own 116-setting Hall database; the rest
// of the search is shared with the 3D space groups, so the family is carried as
// a compile-time parameter and every family-dependent branch is `if constexpr`.
enum class GroupFamily { space, layer };

// Number of Hall settings per family. Static-asserted against the generated
// tables in data/spg_database.hpp.
inline constexpr int kSpaceHallSettings = 530;
inline constexpr int kLayerHallSettings = 116;

[[nodiscard]] constexpr int hall_settings(GroupFamily family) noexcept {
  return family == GroupFamily::layer ? kLayerHallSettings : kSpaceHallSettings;
}

// Bridge to the database's signed-Hall convention (3D settings 1..530, layer
// settings -1..-116), which the catalogs still use. Phase 3 replaces this with
// a validated HallNumber{family, index} key and the convention disappears.
[[nodiscard]] constexpr int signed_hall(GroupFamily family,
                                        int index) noexcept {
  return family == GroupFamily::layer ? -index : index;
}

[[nodiscard]] constexpr GroupFamily family_of_hall(int signed_hall) noexcept {
  return signed_hall < 0 ? GroupFamily::layer : GroupFamily::space;
}

} // namespace cppcrystal
