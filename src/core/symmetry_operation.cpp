#include <seitz/core/symmetry_operation.hpp>

#include "math/integer_matrix.hpp"

namespace seitz {

std::optional<SymmetryOperation> SymmetryOperation::inverse() const noexcept {
  auto const rinv = math::integer_inverse(rotation);
  if (!rinv) {
    return std::nullopt;
  }
  return SymmetryOperation{
      .rotation = *rinv, .translation = -(rinv->cast<double>() * translation)};
}

} // namespace seitz
