#pragma once

#include <cppcrystal/core/types.hpp>

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

// Entry-wise tolerance test against zero: every |entry| strictly below tol.
[[nodiscard]] bool approx_zero(MatrixExpr auto const &a, double tol) noexcept {
  return a.cwiseAbs().maxCoeff() < tol;
}

// Entry-wise tolerance comparison: every |a - b| entry strictly below tol.
// The single definition of "equal within tolerance" for vectors and matrices.
//
// Strictly below, not at: a handful of spglib-parity sites deliberately accept
// exactly-at-tolerance inputs and so spell their own `<= tol` (same_operation
// in core/symmetry_operation.hpp, matrices_close in spacegroup/spacegroup.cpp,
// the lattice check in refine/operations.cpp, and fixes() in
// group/point_group.cpp, whose `none_of(... > kTol)` is the same inclusive
// test turned inside out). Those are not oversights -- do not fold them in
// here without an oracle run, because the two differ at exactly one input
// value and that value is reachable.
[[nodiscard]] bool approx_equal(MatrixExpr auto const &a,
                                MatrixExpr auto const &b, double tol) noexcept {
  return approx_zero(a - b, tol);
}

// An angle tolerance in degrees, or std::nullopt to derive an effective value
// from symprec.
using AngleTolerance = std::optional<double>;

// Bundle of tolerances threaded through the symmetry-finding pipeline: the one
// parameter every stage takes, in place of the (symprec, angle) pair.
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

} // namespace cppcrystal
