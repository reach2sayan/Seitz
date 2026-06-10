#pragma once

#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index_container.hpp>

#include <string_view>

// Access to the built-in space-group database (530 Hall settings). The raw data
// tables are generated into data/spacegroup_tables.hpp by
// tools/transcribe_spg_database.py; here we decode them once into a multi-keyed
// catalog (metadata) plus a flat_map (operations). Port of the spgdb_*
// functions in spg_database.c. Layer groups are out of scope.
namespace spglib::data {

// Bravais centering (spglib Centering enum).
enum class Centering {
  error,
  primitive,
  body,
  face,
  a_face,
  b_face,
  c_face,
  base,
  r_center,
};

struct SpacegroupType {
  int hall_number = 0; // 1..530 (the catalog's unique key)
  int number = 0;      // international space-group number 1..230
  std::string_view schoenflies;
  std::string_view hall_symbol;
  std::string_view international;
  std::string_view international_full;
  std::string_view international_short;
  std::string_view choice;
  Centering centering = Centering::error;
  int pointgroup_number = 0; // 1..32
};

// Index tags for the catalog.
struct ByHall {};
struct ByNumber {};
struct ByPointgroup {};

namespace mi = boost::multi_index;

// Multi-keyed catalog of the 530 Hall settings, mirroring Forcesmith's
// ElementCatalog idiom. Look up a unique setting by Hall number, or enumerate
// every setting of a given space-group number / point group.
using SpacegroupCatalog = mi::multi_index_container<
    SpacegroupType,
    mi::indexed_by<mi::hashed_unique<mi::tag<ByHall>,
                                     mi::member<SpacegroupType, int,
                                                &SpacegroupType::hall_number>>,
                   mi::ordered_non_unique<mi::tag<ByNumber>,
                                          mi::member<SpacegroupType, int,
                                                     &SpacegroupType::number>>,
                   mi::ordered_non_unique<
                       mi::tag<ByPointgroup>,
                       mi::member<SpacegroupType, int,
                                  &SpacegroupType::pointgroup_number>>>>;

struct Database {
  SpacegroupCatalog catalog; // metadata, multi-keyed
  boost::container::flat_map<int, SymmetryOperations> operations; // hall -> ops
};

inline constexpr int kNumHallNumbers = 530;

// The lazily-built, immutable space-group database.
[[nodiscard]] Database const &database();

// Convenience lookups by Hall number (1..530); empty / number-0 if out of
// range.
[[nodiscard]] SymmetryOperations operations_from_database(int hall_number);
[[nodiscard]] SpacegroupType spacegroup_type(int hall_number);

} // namespace spglib::data
