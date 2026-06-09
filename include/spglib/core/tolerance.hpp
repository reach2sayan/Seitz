#pragma once

#include <optional>

namespace spglib {

// Default Cartesian distance tolerance for the symmetry search (spglib
// `symprec`).
inline constexpr double kDefaultSymprec = 1e-5;

// Threshold used when folding fractional coordinates into the unit cell
// (spglib's ZERO_PREC in mathfunc.c). Keeps values a hair below 0 near 0
// instead of wrapping them up to ~1.
inline constexpr double kZeroPrec = 1e-10;

// An angle tolerance in degrees, or std::nullopt to derive an effective value
// from symprec (spglib expresses the latter as a negative angle tolerance).
using AngleTolerance = std::optional<double>;

// Bundle of tolerances threaded through the symmetry-finding pipeline.
struct Tolerance {
  double symprec = kDefaultSymprec;
  AngleTolerance angle_tolerance;
};

} // namespace spglib
