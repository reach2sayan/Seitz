#pragma once

#include <Eigen/Core>

#include <optional>

namespace cppcrystal {

// Default Cartesian distance tolerance for the symmetry search.
inline constexpr double kDefaultSymprec = 1e-5;

// Threshold used when folding fractional coordinates into the unit cell. Keeps
// values a hair below 0 near 0 instead of wrapping them up to ~1.
inline constexpr double kZeroPrec = 1e-10;

// True when |lhs| is longer than |rhs| by more than the near-tie tolerance,
// given their squared norms. Order is preserved on near-equal lengths.
[[nodiscard]] constexpr bool sqnorm_longer(double lhs_sqnorm,
                                           double rhs_sqnorm) {
  return lhs_sqnorm > rhs_sqnorm + kZeroPrec;
}

// Entry-wise tolerance comparison: every |a - b| entry strictly below tol.
// The single definition of "equal within tolerance" for vectors and matrices.
template <typename A, typename B>
[[nodiscard]] bool approx_equal(Eigen::MatrixBase<A> const &a,
                                Eigen::MatrixBase<B> const &b,
                                double tol) noexcept {
  return ((a - b).cwiseAbs().maxCoeff() < tol);
}

// An angle tolerance in degrees, or std::nullopt to derive an effective value
// from symprec.
using AngleTolerance = std::optional<double>;

// Bundle of tolerances threaded through the symmetry-finding pipeline.
struct Tolerance {
  double symprec = kDefaultSymprec;
  AngleTolerance angle_tolerance;
};

} // namespace cppcrystal
