#pragma once

#include "math/integer_matrix.hpp"
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/data/spg_database.hpp> // data::Centering

#include <array>
#include <cstddef>
#include <span>

// Centering change-of-basis data, the single source of truth shared by the
// space-group search (spacegroup.cpp / hall_symbol.cpp) and cell
// standardization (standardize_cell). For a centering C:
//   - centering_matrix(C)      M     : primitive -> conventional basis
//   (integer)
//   - centering_matrix_inv(C)  M^-1  : conventional -> primitive basis
//   - centering_shifts(C)            : the non-trivial centering translations
//                                      (multiplicity - 1 of them; the zero
//                                      shift is implicit)
// PRIMITIVE (and any unused entry) maps to the identity / an empty shift list.
//
// M is the only transcribed table. M^-1 is *derived* from it at compile time as
// adjugate(M) / det(M), which is exact for every centering, so there is no
// second table to keep in sync and no M . M^-1 = I transcription check to run --
// the identity holds by construction (a static_assert below pins that claim).
// The table is a constexpr std::array proxy, row-major, because Eigen types are
// literal for construction and access but not for arithmetic; Eigen
// materialises once at first use.
namespace cppcrystal {

inline constexpr std::size_t kCenteringCount = 9; // Centering::error..r_center

namespace detail {

[[nodiscard]] consteval std::size_t centering_index(data::Centering c) {
  return static_cast<std::size_t>(c);
}

// Row-major 3x3, indexed by Centering. Unset entries are the identity.
inline constexpr auto kCenteringMatrix = [] {
  std::array<math::Mat3Rows, kCenteringCount> t{};
  for (auto &m : t) {
    m = {1, 0, 0, 0, 1, 0, 0, 0, 1};
  }
  t[centering_index(data::Centering::body)] = {0, 1, 1, 1, 0, 1, 1, 1, 0};
  t[centering_index(data::Centering::face)] = {-1, 1, 1, 1, -1, 1, 1, 1, -1};
  t[centering_index(data::Centering::a_face)] = {1, 0, 0, 0, 1, 1, 0, -1, 1};
  t[centering_index(data::Centering::b_face)] = {1, 0, 1, 0, 1, 0, -1, 0, 1};
  t[centering_index(data::Centering::c_face)] = {1, -1, 0, 1, 1, 0, 0, 0, 1};
  t[centering_index(data::Centering::r_center)] = {1, 0, 1, -1, 1, 1, 0, -1, 1};
  return t;
}();

// M^-1 = adjugate(M) / det(M), exact for every centering: det is 1, 2, 3 or 4
// and the adjugate entries are divisible into it in binary floating point
// without rounding. Kept unreduced -- reducing by gcd would change num and den
// by the same factor and land on the identical double.
inline constexpr auto kCenteringMatrixInvNum = [] {
  std::array<math::Mat3Rows, kCenteringCount> t{};
  for (std::size_t i = 0; i < kCenteringCount; ++i) {
    t[i] = math::adjugate(kCenteringMatrix[i]);
  }
  return t;
}();

inline constexpr auto kCenteringMatrixInvDen = [] {
  std::array<int, kCenteringCount> d{};
  for (std::size_t i = 0; i < kCenteringCount; ++i) {
    d[i] = math::determinant(kCenteringMatrix[i]);
  }
  return d;
}();

// The derivation is only worth trusting if M . adjugate(M) == det(M) . I, so
// assert exactly that -- on the derived values, not on transcribed ones.
[[nodiscard]] consteval bool inverses_are_exact() {
  for (std::size_t c = 0; c < kCenteringCount; ++c) {
    auto const &m = kCenteringMatrix[c];
    auto const &n = kCenteringMatrixInvNum[c];
    int const den = kCenteringMatrixInvDen[c];
    if (den == 0) {
      return false; // singular M: no inverse to derive
    }
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        int sum = 0;
        for (std::size_t k = 0; k < 3; ++k) {
          sum += m[i * 3 + k] * n[k * 3 + j];
        }
        if (sum != (i == j ? den : 0)) {
          return false;
        }
      }
    }
  }
  return true;
}
static_assert(inverses_are_exact());

} // namespace detail

[[nodiscard]] Matrix3i const &centering_matrix(data::Centering c);
[[nodiscard]] Matrix3d const &centering_matrix_inv(data::Centering c);
[[nodiscard]] std::span<Vector3d const> centering_shifts(data::Centering c);

} // namespace cppcrystal
