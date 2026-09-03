#include "core/centering.hpp"

#include <array>
#include <cstddef>
#include <span>

namespace cppcrystal {

using data::Centering;

namespace {

// The proxies materialise into Eigen once, on first use: the tables are the
// source of truth, these are just the typed views the pipeline works with.
template <class Matrix, class Proxy>
[[nodiscard]] Matrix materialise(Proxy const &rows, auto scale) {
  Matrix m;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      m(i, j) = scale(rows[static_cast<std::size_t>(i * 3 + j)]);
    }
  }
  return m;
}

} // namespace

Matrix3i const &centering_matrix(Centering c) {
  static auto const table = [] {
    std::array<Matrix3i, kCenteringCount> t;
    for (std::size_t i = 0; i < kCenteringCount; ++i) {
      t[i] = materialise<Matrix3i>(detail::kCenteringMatrix[i],
                                   [](int v) { return v; });
    }
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

Matrix3d const &centering_matrix_inv(Centering c) {
  static auto const table = [] {
    std::array<Matrix3d, kCenteringCount> t;
    for (std::size_t i = 0; i < kCenteringCount; ++i) {
      double const den =
          static_cast<double>(detail::kCenteringMatrixInvDen[i]);
      t[i] = materialise<Matrix3d>(
          detail::kCenteringMatrixInvNum[i],
          [den](int v) { return static_cast<double>(v) / den; });
    }
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

} // namespace cppcrystal
