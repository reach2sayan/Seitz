#pragma once

#include <spglib/core/types.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>

#include <optional>
#include <vector>

namespace spglib {

// A space-group symmetry operation acting on fractional coordinates:
//   x -> rotation . x + translation.
struct SymmetryOperation {
  Matrix3i rotation{Matrix3i::Identity()};
  Vector3d translation{Vector3d::Zero()};

  [[nodiscard]] Vector3d apply(Vector3d const &x) const noexcept {
    return rotation.cast<double>() * x + translation;
  }

  // Group composition: (a * b).apply(x) == a.apply(b.apply(x)).
  [[nodiscard]] SymmetryOperation
  operator*(SymmetryOperation const &rhs) const noexcept {
    return {rotation * rhs.rotation,
            rotation.cast<double>() * rhs.translation + translation};
  }

  // Inverse operation; std::nullopt when the rotation is not unimodular.
  [[nodiscard]] std::optional<SymmetryOperation> inverse() const noexcept {
    auto const rinv = math::integer_inverse(rotation);
    if (!rinv) {
      return std::nullopt;
    }
    return SymmetryOperation{*rinv, -(rinv->cast<double>() * translation)};
  }

  [[nodiscard]] bool is_identity_rotation() const noexcept {
    return rotation == Matrix3i::Identity();
  }
};

using SymmetryOperations = std::vector<SymmetryOperation>;
[[nodiscard]] inline bool same_operation(SymmetryOperation const &a,
                                         SymmetryOperation const &b,
                                         double symprec) noexcept {
  if (a.rotation != b.rotation) {
    return false;
  }
  Vector3d const d = math::nearest_offset(Vector3d(a.translation - b.translation));
  return d.cwiseAbs().maxCoeff() <= symprec;
}

} // namespace spglib
