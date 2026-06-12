#include <spglib/data/rod_database.hpp>

#include <spglib/data/rod_group_tables.hpp>

namespace spglib::data {

namespace {

[[nodiscard]] bool in_range(int rod_number) noexcept {
  return rod_number >= 1 && rod_number <= kNumRodGroups;
}

} // namespace

int num_rod_groups() noexcept { return kNumRodGroups; }

int rod_periodic_axis() noexcept { return kRodPeriodicAxis; }

std::string_view rod_symbol(int rod_number) noexcept {
  if (!in_range(rod_number)) {
    return {};
  }
  return kRodGroupSymbols[static_cast<std::size_t>(rod_number - 1)];
}

SymmetryOperations rod_operations_from_database(int rod_number) {
  if (!in_range(rod_number)) {
    return {};
  }
  int const begin = kRodOperationOffset[static_cast<std::size_t>(rod_number - 1)];
  int const end = kRodOperationOffset[static_cast<std::size_t>(rod_number)];

  SymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(end - begin));
  for (int i = begin; i < end; ++i) {
    RodOperation const &raw = kRodOperations[static_cast<std::size_t>(i)];

    Matrix3i rotation;
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        rotation(r, c) = raw.rotation[static_cast<std::size_t>(r * 3 + c)];
      }
    }
    Vector3d translation;
    for (int a = 0; a < 3; ++a) {
      translation[a] = static_cast<double>(raw.translation[static_cast<std::size_t>(a)]) /
                       static_cast<double>(kRodTranslationDenominator);
    }
    ops.push_back(SymmetryOperation{rotation, translation});
  }
  return ops;
}

} // namespace spglib::data
