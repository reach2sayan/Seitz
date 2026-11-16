#include "data/rod_database.hpp"
#include <seitz/core/operation_set.hpp>

#include "data/rod_group_tables.hpp"
#include "math/integer_matrix.hpp"

#include <ranges>
#include <span>

namespace seitz::data {

Operations rod_operations_from_database(int rod_number) {
  auto const begin = static_cast<std::size_t>(
      kRodOperationOffset[static_cast<std::size_t>(rod_number - 1)]);
  auto const end = static_cast<std::size_t>(
      kRodOperationOffset[static_cast<std::size_t>(rod_number)]);
  return Operations{
      std::from_range,
      std::span{kRodOperations}.subspan(begin, end - begin) |
          std::views::transform([](RodOperation const &raw) {
            return SymmetryOperation{
                .rotation = math::as_matrix(raw.rotation),
                .translation =
                    Eigen::Map<Eigen::Vector3i const>(raw.translation.data())
                        .cast<double>() /
                    static_cast<double>(kRodTranslationDenominator)};
          })};
}

} // namespace seitz::data
