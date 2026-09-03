#include <cppcrystal/core/symmetry_operation.hpp>

#include "math/integer_matrix.hpp"

namespace cppcrystal {

std::optional<SymmetryOperation> SymmetryOperation::inverse() const noexcept {
  auto const rinv = math::integer_inverse(rotation);
  if (!rinv) {
    return std::nullopt;
  }
  return SymmetryOperation{*rinv, -(rinv->cast<double>() * translation)};
}

} // namespace cppcrystal
