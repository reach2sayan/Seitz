#include "core/centering.hpp"

#include <array>
#include <span>
#include <cstddef>

namespace cppcrystal {

using data::Centering;

namespace {
constexpr std::size_t kCenteringCount = 9; // Centering::error .. r_center
} // namespace

Matrix3i const &centering_matrix(Centering c) {
  static auto const table = [] {
    std::array<Matrix3i, kCenteringCount> t;
    t.fill(Matrix3i::Identity());
    auto at = [&](Centering cc) -> Matrix3i & {
      return t[static_cast<std::size_t>(cc)];
    };
    at(Centering::body) << 0, 1, 1, 1, 0, 1, 1, 1, 0;
    at(Centering::face) << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    at(Centering::a_face) << 1, 0, 0, 0, 1, 1, 0, -1, 1;
    at(Centering::b_face) << 1, 0, 1, 0, 1, 0, -1, 0, 1;
    at(Centering::c_face) << 1, -1, 0, 1, 1, 0, 0, 0, 1;
    at(Centering::r_center) << 1, 0, 1, -1, 1, 1, 0, -1, 1;
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

Matrix3d const &centering_matrix_inv(Centering c) {
  static auto const table = [] {
    std::array<Matrix3d, kCenteringCount> t;
    t.fill(Matrix3d::Identity());
    auto at = [&](Centering cc) -> Matrix3d & {
      return t[static_cast<std::size_t>(cc)];
    };
    at(Centering::body) << -0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5;
    at(Centering::face) << 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0;
    at(Centering::a_face) << 1.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5;
    at(Centering::b_face) << 0.5, 0.0, -0.5, 0.0, 1.0, 0.0, 0.5, 0.0, 0.5;
    at(Centering::c_face) << 0.5, 0.5, 0.0, -0.5, 0.5, 0.0, 0.0, 0.0, 1.0;
    at(Centering::r_center) << 2. / 3, -1. / 3, -1. / 3, 1. / 3, 1. / 3, -2. / 3,
        1. / 3, 1. / 3, 1. / 3;
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

std::span<Vector3d const> centering_shifts(Centering centering) {
  static std::array<Vector3d, 1> const a{{{0.0, 0.5, 0.5}}};
  static std::array<Vector3d, 1> const b{{{0.5, 0.0, 0.5}}};
  static std::array<Vector3d, 1> const c{{{0.5, 0.5, 0.0}}};
  static std::array<Vector3d, 1> const body{{{0.5, 0.5, 0.5}}};
  static std::array<Vector3d, 2> const rhombohedral{
      {{2. / 3, 1. / 3, 1. / 3}, {1. / 3, 2. / 3, 2. / 3}}};
  static std::array<Vector3d, 3> const face{
      {{0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}}};
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
