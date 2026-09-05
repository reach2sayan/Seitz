#pragma once

#include <seitz/core/types.hpp>

#include <optional>

#pragma GCC visibility push(default)

namespace seitz {

// Default Cartesian distance tolerance for the symmetry search.
inline constexpr double kDefaultSymprec = 1e-5;

// Fold threshold: keeps values a hair below 0 near 0, not wrapped to ~1.
inline constexpr double kZeroPrec = 1e-10;

// |lhs| > |rhs| by more than the near-tie tolerance, from squared norms;
// near-equal lengths keep their order.
[[nodiscard]] constexpr bool sqnorm_longer(double lhs_sqnorm,
                                           double rhs_sqnorm) {
  return lhs_sqnorm > rhs_sqnorm + kZeroPrec;
}

// Entry-wise tolerance test against zero: every |entry| strictly below tol.
[[nodiscard]] bool approx_zero(MatrixExpr auto const &a, double tol) noexcept {
  return a.cwiseAbs().maxCoeff() < tol;
}

// Entry-wise tolerance comparison: every |a - b| entry strictly below tol.
[[nodiscard]] bool approx_equal(MatrixExpr auto const &a,
                                MatrixExpr auto const &b, double tol) noexcept {
  return approx_zero(a - b, tol);
}

// An angle tolerance in degrees, or std::nullopt to derive an effective value
// from symprec.
using AngleTolerance = std::optional<double>;

// Bundle of tolerances threaded through the symmetry-finding pipeline:
struct Tolerance {
  double symprec = kDefaultSymprec;
  AngleTolerance angle_tolerance;
};

// The magnetic search adds a moment-comparison tolerance, which defaults to
// symprec when unset.
struct MagneticTolerance : Tolerance {
  std::optional<double> moment;
  [[nodiscard]] double moment_or_symprec() const noexcept {
    return moment.value_or(symprec);
  }
};

} // namespace seitz

#pragma GCC visibility pop
