#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/data/spg_database.hpp>

#include <cppcrystal/data/detail/lookup.hpp>
#include <cppcrystal/data/detail/packed_decode.hpp>
#include <cppcrystal/data/spacegroup_operation_tables.hpp>

#include <array>

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

Operations const &operations_from_database(HallNumber hall) {
  // Materialised once per family. The generated index tables keep their 1-based
  // layout, so a HallNumber's index addresses them directly.
  static auto const ops =
      detail::build_operation_table(kSymmetryOperationIndex, kDecodedOps);
  static auto const layer_ops =
      detail::build_operation_table(kLayerSymmetryOperationIndex, kDecodedOps);
  auto const i = static_cast<std::size_t>(hall.index());
  return hall.family() == GroupFamily::layer ? layer_ops[i] : ops[i];
}

} // namespace cppcrystal::data
