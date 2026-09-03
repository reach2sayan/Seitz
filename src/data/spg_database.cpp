#include <cppcrystal/data/spg_database.hpp>

#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/detail/packed_decode.hpp>
#include <cppcrystal/data/spacegroup_operation_tables.hpp>

#include <array>
#include <vector>

namespace cppcrystal::data {

namespace {

// The whole packed table, unpacked once at compile time so no base-3/12
// decoding happens at runtime.
constexpr auto kDecodedOps = [] {
  std::array<detail::DecodedOp, kSymmetryOperations.size()> out{};
  for (std::size_t k = 0; k < kSymmetryOperations.size(); ++k) {
    out[k] = detail::decode_packed(kSymmetryOperations[k]);
  }
  return out;
}();

} // namespace

SymmetryOperations const &operations_from_database(int hall_number) {
  // Materialised once per family, indexed by Hall number (layer settings by
  // -hall); entry 0 is the empty fallback for out-of-range queries.
  static auto const ops =
      detail::build_operation_table(kSymmetryOperationIndex, kDecodedOps);
  static auto const layer_ops =
      detail::build_operation_table(kLayerSymmetryOperationIndex, kDecodedOps);
  return detail::hall_indexed(ops, layer_ops, hall_number);
}

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
