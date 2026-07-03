#include <cppcrystal/core/centering.hpp>

#include <array>
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

std::vector<Vector3d> centering_shifts(Centering centering) {
  switch (centering) {
  case Centering::a_face:
    return {{0.0, 0.5, 0.5}};
  case Centering::b_face:
    return {{0.5, 0.0, 0.5}};
  case Centering::c_face:
    return {{0.5, 0.5, 0.0}};
  case Centering::body:
    return {{0.5, 0.5, 0.5}};
  case Centering::r_center:
    return {{2. / 3, 1. / 3, 1. / 3}, {1. / 3, 2. / 3, 2. / 3}};
  case Centering::face:
    return {{0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
  default:
    return {};
  }
}

} // namespace cppcrystal
