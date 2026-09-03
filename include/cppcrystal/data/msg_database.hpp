#pragma once

#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/magnetic_spacegroup_metadata_tables.hpp>
#include <cppcrystal/data/spg_database.hpp> // kNumHallNumbers

#include <array>
#include <cstddef>
#include <optional>
#include <string_view>

// Access to the built-in magnetic space-group database (1651 UNI numbers). The
// raw data tables are generated into two headers:
// magnetic_spacegroup_metadata_tables.hpp (catalog metadata + the small
// Hall<->UNI mapping tables, included here to back the constexpr catalog) and
// magnetic_spacegroup_operation_tables.hpp (the encoded operations +
// alternative transformations — a compile-time-only compaction, included only
// by the .cpp). As with spg_database, metadata is decoded once at compile time
// into a constexpr catalog, while the symmetry operations (which carry Eigen,
// not a literal type) are decoded on demand.
namespace cppcrystal::data {

inline constexpr int kNumUniNumbers = 1651;

// Magnetic space-group type metadata.
struct MagneticSpacegroupType {
  int uni_number = 0;          // 1..1651 (the catalog's unique key)
  int litvin_number = 0;       // Litvin's sequential number
  std::string_view bns_number; // Belov-Neronova-Smirnova symbol
  std::string_view og_number;  // Opechowski-Guccione symbol
  int number = 0;              // family-space-group international number 1..230
  int type = 0;                // construction type 1..4 (type-I..IV)
};

struct MagneticSpacegroupCatalog {
  std::array<MagneticSpacegroupType, kNumUniNumbers + 1> by_uni{};

  [[nodiscard]] constexpr MagneticSpacegroupType const &
  at(int uni) const noexcept {
    return detail::at_or_sentinel(by_uni, uni);
  }
};

inline constexpr MagneticSpacegroupCatalog kMagneticCatalog = [] {
  MagneticSpacegroupCatalog c{};
  for (std::size_t u = 0; u < kMagneticSpacegroupTypes.size(); ++u) {
    auto const &r = kMagneticSpacegroupTypes[u];
    c.by_uni[u] =
        MagneticSpacegroupType{r.uni_number, r.litvin_number, r.bns_number,
                               r.og_number,  r.number,        r.type};
  }
  return c;
}();

[[nodiscard]] constexpr MagneticSpacegroupType
magnetic_spacegroup_type(int uni_number) noexcept {
  return kMagneticCatalog.at(uni_number);
}

[[nodiscard]] constexpr std::optional<std::pair<int, int>>
uni_candidates(int hall_number) noexcept {
  if (hall_number < 1 || hall_number > kNumHallNumbers) {
    return std::nullopt;
  }
  auto const &m = kMagneticHallMapping[static_cast<std::size_t>(hall_number)];
  return std::pair<int, int>{m[0], m[1]};
}

// The operations of a UNI number in the given Hall setting (0 = its first
// setting), materialised once per setting; empty when (uni, hall) is not a
// valid pairing.
[[nodiscard]] MagneticSymmetryOperations const &
magnetic_operations_from_database(int uni_number, int hall_number = 0);

// The alternative standardized-setting transformations of a UNI number in the
// given Hall setting, identity first; empty when (uni, hall) is not valid.
[[nodiscard]] SymmetryOperations const &
magnetic_std_transformations(int uni_number, int hall_number = 0);

} // namespace cppcrystal::data
