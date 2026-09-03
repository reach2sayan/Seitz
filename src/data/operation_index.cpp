#include "data/operation_index.hpp"

#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/spacegroup_operation_tables.hpp>

#include <cstddef>
#include <vector>

namespace cppcrystal::data {

RotationMultimap<int> const &operations_by_rotation(int hall_number) {
  // Entry 0 (the empty fallback) indexes the empty operation list.
  auto const index_family = [](int sign, std::size_t size) {
    std::vector<RotationMultimap<int>> table(size);
    for (std::size_t key = 1; key < size; ++key) {
      table[key] = index_by_rotation(
          operations_from_database(sign * static_cast<int>(key)),
          &SymmetryOperation::rotation);
    }
    return table;
  };
  static auto const index = index_family(1, kSymmetryOperationIndex.size());
  static auto const layer_index =
      index_family(-1, kLayerSymmetryOperationIndex.size());
  return detail::hall_indexed(index, layer_index, hall_number);
}

} // namespace cppcrystal::data
