#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/data/catalog.hpp>
#include <seitz/data/detail/lookup.hpp>
#include <seitz/data/magnetic_spacegroup_metadata_tables.hpp>
#include <seitz/data/spg_database.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <string_view>

// The built-in magnetic space-group database (1651 UNI numbers). Raw tables are
// generated into magnetic_spacegroup_metadata_tables.hpp (metadata + the
// Hall<->UNI mappings, backing the constexpr catalog here) and
// magnetic_spacegroup_operation_tables.hpp (encoded operations + alternative
// transformations, included only by the .cpp). As in spg_database, metadata
// decodes at compile time, operations on demand.
#pragma GCC visibility push(default)

namespace seitz::data {

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

// The magnetic family policy: 1651 UNI numbers, keyed by UniNumber.
struct MagneticFamily {
  using Row = MagneticSpacegroupType;
  using Key = UniNumber;
  static constexpr std::size_t count = static_cast<std::size_t>(kUniNumbers);
  [[nodiscard]] static constexpr int index_of(UniNumber key) noexcept {
    return key.value();
  }

  // The generated table is 1-based with a dummy row 0; the catalog is not.
  [[nodiscard]] static constexpr std::array<Row, count> decode() {
    std::array<Row, count> rows{};
    for (std::size_t i = 0; i < count; ++i) {
      auto const &r = kMagneticSpacegroupTypes[i + 1];
      rows[i] =
          MagneticSpacegroupType{r.uni_number, r.litvin_number, r.bns_number,
                                 r.og_number,  r.number,        r.type};
    }
    return rows;
  }
};

[[nodiscard]] constexpr MagneticSpacegroupType const &
magnetic_spacegroup_type(UniNumber uni) noexcept {
  return kCatalog<MagneticFamily>[uni];
}

// The [first, last] UNI numbers a 3D Hall setting can carry.
[[nodiscard]] constexpr std::pair<UniNumber, UniNumber>
uni_candidates(HallNumber hall) noexcept {
  auto const &m = kMagneticHallMapping[static_cast<std::size_t>(hall.index())];
  return {*UniNumber::of(m[0]), *UniNumber::of(m[1])};
}

// The operations of a UNI number in the given Hall setting (unset = its first
// setting), materialised once per setting; empty when (uni, hall) is not a
// valid pairing.
[[nodiscard]] MagneticOperations const &magnetic_operations_from_database(
    UniNumber uni, std::optional<HallNumber> hall = std::nullopt);

// The alternative standardized-setting transformations of a UNI number in the
// given Hall setting, identity first; empty when (uni, hall) is not valid.
[[nodiscard]] Operations const &
magnetic_std_transformations(UniNumber uni,
                             std::optional<HallNumber> hall = std::nullopt);

} // namespace seitz::data

#pragma GCC visibility pop
