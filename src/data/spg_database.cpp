#include <spglib/data/spg_database.hpp>

#include <spglib/data/spacegroup_operation_tables.hpp>

#include <array>
#include <cstdint>

namespace spglib::data {

namespace {

// Decoded form of one packed operation: rotation as 9 row-major elements in
// {-1,0,1}, translation as 3 base-12 numerators (value n/12). This is spglib's
// spgdb_decode_symmetry encoding (rotation base-3, translation base-12),
// unpacked once at compile time so no decoding happens at runtime.
struct DecodedOp {
  std::array<std::int8_t, 9> rot;
  std::array<std::int8_t, 3> trans_num;
};

constexpr auto kDecodedOps = [] {
  std::array<DecodedOp, kSymmetryOperations.size()> out{};
  for (std::size_t k = 0; k < kSymmetryOperations.size(); ++k) {
    int const encoded = kSymmetryOperations[k];
    int const r = encoded % 19683; // 3^9
    int digit = 6561;              // 3^8
    for (std::size_t e = 0; e < 9; ++e) {
      out[k].rot[e] = static_cast<std::int8_t>((r % (digit * 3)) / digit - 1);
      digit /= 3;
    }
    int const t = encoded / 19683;
    digit = 144; // 12^2
    for (std::size_t e = 0; e < 3; ++e) {
      out[k].trans_num[e] =
          static_cast<std::int8_t>((t % (digit * 12)) / digit);
      digit /= 12;
    }
  }
  return out;
}();

[[nodiscard]] SymmetryOperation make_operation(DecodedOp const &d) noexcept {
  Matrix3i rot;
  rot << d.rot[0], d.rot[1], d.rot[2], d.rot[3], d.rot[4], d.rot[5], d.rot[6],
      d.rot[7], d.rot[8];

  Vector3d const trans{d.trans_num[0] / 12.0, d.trans_num[1] / 12.0,
                       d.trans_num[2] / 12.0};

  return {rot, trans};
}

// Materialise every Hall setting's operations once, indexed by Hall number
// (1..530). Index 0 is the empty sentinel for out-of-range queries.
[[nodiscard]] std::array<SymmetryOperations, kNumHallNumbers + 1>
build_operations() {
  std::array<SymmetryOperations, kNumHallNumbers + 1> ops;
  for (int hall = 1; hall <= kNumHallNumbers; ++hall) {
    auto const &idx = kSymmetryOperationIndex[static_cast<std::size_t>(hall)];
    SymmetryOperations v;
    v.reserve(static_cast<std::size_t>(idx[0]));
    for (int i = 0; i < idx[0]; ++i) {
      v.push_back(
          make_operation(kDecodedOps[static_cast<std::size_t>(idx[1] + i)]));
    }
    ops[static_cast<std::size_t>(hall)] = std::move(v);
  }
  return ops;
}

} // namespace

SymmetryOperations const &operations_from_database(int hall_number) {
  static std::array<SymmetryOperations, kNumHallNumbers + 1> const ops =
      build_operations();
  bool const in_range = hall_number >= 1 && hall_number <= kNumHallNumbers;
  return ops[static_cast<std::size_t>(in_range ? hall_number : 0)];
}

} // namespace spglib::data
