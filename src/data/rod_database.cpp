#include <cppcrystal/core/operation_set.hpp>
#include "data/rod_database.hpp"

#include <cppcrystal/data/rod_group_tables.hpp>

namespace cppcrystal::data {

Operations rod_operations_from_database(int rod_number) {
  if (!rod_number_in_range(rod_number)) {
    return {};
  }
  int const begin = kRodOperationOffset[static_cast<std::size_t>(rod_number - 1)];
  int const end = kRodOperationOffset[static_cast<std::size_t>(rod_number)];

  std::vector<SymmetryOperation> ops;
  ops.reserve(static_cast<std::size_t>(end - begin));
  for (int i = begin; i < end; ++i) {
    RodOperation const &raw = kRodOperations[static_cast<std::size_t>(i)];

    Matrix3i rotation =
        Eigen::Map<Eigen::Matrix<int, 3, 3, Eigen::RowMajor> const>(
            raw.rotation.data());
    Vector3d translation =
        Eigen::Map<Eigen::Vector3i const>(raw.translation.data())
            .cast<double>() /
        static_cast<double>(kRodTranslationDenominator);
    ops.push_back({.rotation = rotation, .translation = translation});
  }
  return Operations{std::move(ops)};
}

} // namespace cppcrystal::data
