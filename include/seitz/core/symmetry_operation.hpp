#pragma once

#include <seitz/core/fractional.hpp>
#include <seitz/core/types.hpp>

#include <concepts>
#include <optional>

#pragma GCC visibility push(default)

namespace seitz {

// Anything shaped like a space-group operation: an integer rotation plus a
// fractional translation (SymmetryOperation, MagneticSymmetryOperation).
// Constrains the generic algorithms that only touch the spatial part.
template <class Op>
concept SpaceGroupOperationLike = requires(Op op) {
  { op.rotation } -> std::convertible_to<Matrix3i>;
  { op.translation } -> std::convertible_to<Vector3d>;
};

// Whether an operation set includes its time-reversal partners: the magnetic
// search's anti-operations, and the reciprocal point group's inversion partner.
enum class TimeReversal { off, on };

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

  // Inverse operation; std::nullopt when the rotation is not unimodular. Out
  // of line: the exact integer inverse it needs is private to src/math.
  [[nodiscard]] std::optional<SymmetryOperation> inverse() const noexcept;

  [[nodiscard]] bool is_identity_rotation() const noexcept {
    return rotation == Matrix3i::Identity();
  }
};

// Change of basis of an operation: (T, 0)(R, t)(T, 0)^-1 = (T R T^-1, T t),
// with the conjugated rotation rounded back to the integer basis. Any extra
// fields of `Op` (e.g. a time-reversal flag) ride along unchanged.
[[nodiscard]] auto conjugated_by(SpaceGroupOperationLike auto op,
                                 Matrix3d const &t,
                                 Matrix3d const &t_inv) noexcept {
  op.rotation =
      math::round_to_int(t * op.rotation.template cast<double>() * t_inv);
  op.translation = t * op.translation;
  return op;
}
[[nodiscard]] inline bool same_operation(SymmetryOperation const &a,
                                         SymmetryOperation const &b,
                                         double symprec) noexcept {
  if (a.rotation != b.rotation) {
    return false;
  }
  Vector3d const d =
      math::nearest_offset(Vector3d(a.translation - b.translation));
  return d.cwiseAbs().maxCoeff() <= symprec;
}

} // namespace seitz

#pragma GCC visibility pop
