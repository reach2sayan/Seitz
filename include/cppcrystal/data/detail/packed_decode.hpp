#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <Eigen/Dense>

#include <array>
#include <cstddef>
#include <cstdint>

// The packed space-group operation encoding shared by the 3D/layer and the
// magnetic databases: rotation as 9 base-3 digits (values {-1,0,1}), then
// translation as 3 base-12 numerators (value n/12). Decoding is constexpr so
// the tables unpack once at compile time; only the Eigen materialisation runs
// at runtime (Eigen matrices cannot live in constexpr globals).
namespace cppcrystal::data::detail {

// Decoded form of one packed operation: rotation as 9 row-major elements,
// translation as 3 base-12 numerators.
struct DecodedOp {
  std::array<std::int8_t, 9> rot;
  std::array<std::int8_t, 3> trans_num;
};

[[nodiscard]] constexpr DecodedOp decode_packed(int encoded) noexcept {
  DecodedOp d{};
  constexpr int kRotBase = 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3; // 3^9
  int const r = encoded % kRotBase;
  int digit = kRotBase / 3; // 3^8
  for (auto &x : d.rot) {
    x = static_cast<std::int8_t>((r % (digit * 3)) / digit - 1);
    digit /= 3;
  }
  int const t = encoded / kRotBase;
  digit = 12 * 12;
  for (auto &x : d.trans_num) {
    x = static_cast<std::int8_t>((t % (digit * 12)) / digit);
    digit /= 12;
  }
  return d;
}

// Materialise the Eigen-valued operation from decoded data.
[[nodiscard]] inline SymmetryOperation
make_operation(DecodedOp const &d) noexcept {
  Matrix3i const rot =
      Eigen::Map<Eigen::Matrix<std::int8_t, 3, 3, Eigen::RowMajor> const>(
          d.rot.data())
          .cast<int>();
  Vector3d const trans =
      Eigen::Map<Eigen::Matrix<std::int8_t, 3, 1> const>(d.trans_num.data())
          .cast<double>() /
      12.0;
  return {rot, trans};
}

// Materialise per-setting operation lists from a {count, offset} index table
// over a decoded-operation array. Entry 0 of the result stays empty — the
// out-of-range fallback the lookup helpers hand back.
template <std::size_t N, std::size_t M>
[[nodiscard]] std::array<SymmetryOperations, N>
build_operation_table(std::array<std::array<int, 2>, N> const &index,
                      std::array<DecodedOp, M> const &decoded) {
  std::array<SymmetryOperations, N> ops;
  for (std::size_t key = 1; key < N; ++key) {
    auto const [count, offset] = index[key];
    SymmetryOperations v;
    v.reserve(static_cast<std::size_t>(count));
    for (int i = 0; i < count; ++i) {
      v.push_back(make_operation(decoded[static_cast<std::size_t>(offset + i)]));
    }
    ops[key] = std::move(v);
  }
  return ops;
}

} // namespace cppcrystal::data::detail
