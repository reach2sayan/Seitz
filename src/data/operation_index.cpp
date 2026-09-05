#include "data/operation_index.hpp"

#include "data/spacegroup_operation_tables.hpp"
#include <seitz/data/detail/lookup.hpp>

#include <cstddef>
#include <vector>

namespace seitz::data {

RotationMultimap<int> const &operations_by_rotation(HallNumber hall) {
  // 1-based like the generated tables; entry 0 is never addressed.
  auto const index_family = [](GroupFamily family, std::size_t size) {
    std::vector<RotationMultimap<int>> table(size);
    for (std::size_t key = 1; key < size; ++key) {
      table[key] = index_by_rotation(operations_from_database(*HallNumber::of(
                                         family, static_cast<int>(key))),
                                     &SymmetryOperation::rotation);
    }
    return table;
  };
  static auto const index =
      index_family(GroupFamily::space, kSymmetryOperationIndex.size());
  static auto const layer_index =
      index_family(GroupFamily::layer, kLayerSymmetryOperationIndex.size());
  auto const i = static_cast<std::size_t>(hall.index());
  return hall.family() == GroupFamily::layer ? layer_index[i] : index[i];
}

} // namespace seitz::data
