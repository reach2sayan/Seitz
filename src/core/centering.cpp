#include "core/centering.hpp"

#include "math/integer_matrix.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>

namespace seitz {

using data::Centering;

// The tables are the source of truth; these are the typed views the pipeline
// works with, materialised into Eigen once on first use.

Matrix3i const &centering_matrix(Centering c) {
  static auto const table = [] {
    std::array<Matrix3i, kCenteringCount> t;
    std::ranges::transform(detail::kCenteringMatrix, t.begin(),
                           [](math::Mat3Rows const &rows) {
                             return Matrix3i(math::as_matrix(rows));
                           });
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

Matrix3d const &centering_matrix_inv(Centering c) {
  static auto const table = [] {
    std::array<Matrix3d, kCenteringCount> t;
    std::ranges::transform(
        detail::kCenteringMatrixInvNum, detail::kCenteringMatrixInvDen,
        t.begin(), [](math::Mat3Rows const &num, int den) {
          return Matrix3d(math::as_matrix(num).cast<double>() /
                          static_cast<double>(den));
        });
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

std::span<Vector3d const> centering_shifts(Centering centering) {
  // Sixths, so every centering translation (halves and thirds) is exact.
  static std::array<Vector3d, 1> const a{{{0.0, 3. / 6, 3. / 6}}};
  static std::array<Vector3d, 1> const b{{{3. / 6, 0.0, 3. / 6}}};
  static std::array<Vector3d, 1> const c{{{3. / 6, 3. / 6, 0.0}}};
  static std::array<Vector3d, 1> const body{{{3. / 6, 3. / 6, 3. / 6}}};
  static std::array<Vector3d, 2> const rhombohedral{
      {{4. / 6, 2. / 6, 2. / 6}, {2. / 6, 4. / 6, 4. / 6}}};
  static std::array<Vector3d, 3> const face{
      {{0.0, 3. / 6, 3. / 6}, {3. / 6, 0.0, 3. / 6}, {3. / 6, 3. / 6, 0.0}}};
  switch (centering) {
  case Centering::a_face:
    return a;
  case Centering::b_face:
    return b;
  case Centering::c_face:
    return c;
  case Centering::body:
    return body;
  case Centering::r_center:
    return rhombohedral;
  case Centering::face:
    return face;
  default:
    return {};
  }
}

} // namespace seitz
