#pragma once

#include <seitz/core/types.hpp>

#include <algorithm>
#include <cmath>

// The cell-parameter <-> basis-matrix conversion, in the standard
// crystallographic orientation: a along x, b in the xy-plane, c completing a
// right-handed set. Shared by the conventional-cell idealiser (refine) and the
// random-lattice generator (generate), which previously transcribed the same
// two expressions independently.
namespace seitz::math {

// sqrt of the i-th diagonal of a metric (Gram) matrix: |a_i|.
[[nodiscard]] inline double metric_length(Matrix3d const &g, int i) {
  return std::sqrt(g(i, i));
}

// The cosine of the angle between basis vectors i and j of a metric, clamped
// to [-1, 1]. The quotient can exceed 1 by an ulp for a near-degenerate cell,
// where acos would return NaN and sqrt(1 - cos^2) the root of a negative;
// clamping only affects inputs that are already outside the domain, so no
// in-range result changes.
[[nodiscard]] inline double metric_cosine(Matrix3d const &g, int i, int j) {
  return std::clamp(g(i, j) / metric_length(g, i) / metric_length(g, j), -1.0,
                    1.0);
}

// The angle between basis vectors i and j of a metric, in radians.
[[nodiscard]] inline double metric_angle(Matrix3d const &g, int i, int j) {
  return std::acos(metric_cosine(g, i, j));
}

// Basis matrix (columns = basis vectors) from cell parameters, angles in
// RADIANS. The result is upper triangular, which is the same thing as saying
// it is the unique Cholesky factor U of the metric with U^T U = G and a
// positive diagonal.
[[nodiscard]] inline Matrix3d lattice_from_parameters(double a, double b,
                                                      double c, double alpha,
                                                      double beta,
                                                      double gamma) {
  double const ca = std::cos(alpha);
  double const cb = std::cos(beta);
  double const cg = std::cos(gamma);
  double const sg = std::sin(gamma);

  // The discriminant is non-negative for any realisable cell; the clamp only
  // catches a rounding-negative value for a near-degenerate one, where the
  // bare sqrt would be NaN.
  double const cz = c *
                    std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg +
                                                2.0 * ca * cb * cg)) /
                    sg;

  Matrix3d m;
  m.col(0) = Vector3d{a, 0.0, 0.0};
  m.col(1) = Vector3d{b * cg, b * sg, 0.0};
  m.col(2) = Vector3d{c * cb, c * (ca - cb * cg) / sg, cz};
  return m;
}

} // namespace seitz::math
