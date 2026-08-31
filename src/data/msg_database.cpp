#include <cppcrystal/data/msg_database.hpp>

#include <cppcrystal/data/detail/packed_decode.hpp>
#include <cppcrystal/data/magnetic_spacegroup_operation_tables.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <ranges>

namespace cppcrystal::data {

namespace {

using detail::DecodedOp;
using detail::decode_packed;
using detail::make_operation;

// 34012224 = 3^9 * 12^3, the base separating the time-reversal flag from the
// packed spatial operation.
constexpr int kTimeReversalBase =
    (3 * 3 * 3 * 3 * 3 * 3 * 3 * 3 * 3) * (12 * 12 * 12);

// A magnetic operation in decoded literal form: the spatial part plus the
// time-reversal flag.
struct DecodedMagneticOp {
  DecodedOp op;
  bool time_reversal;
};

// Every magnetic operation, decoded once at compile time so no base-3/12
// arithmetic runs at runtime.
constexpr auto kDecodedMagneticOps = [] {
  std::array<DecodedMagneticOp, kMagneticOperations.size()> out{};
  for (std::size_t k = 0; k < kMagneticOperations.size(); ++k) {
    int const enc = kMagneticOperations[k];
    out[k] = {decode_packed(enc % kTimeReversalBase),
              enc / kTimeReversalBase != 0};
  }
  return out;
}();

// Every alternative standardized-setting transformation, decoded once at
// compile time. std::nullopt is the per-row terminator (the source table
// zero-terminates each [uni][offset] row); using optional rather than a magic 0
// keeps the no-sentinel house rule.
constexpr auto kDecodedAltTransformations = [] {
  std::array<std::array<std::array<std::optional<DecodedOp>, 7>, 18>,
             kAlternativeTransformations.size()>
      out{};
  for (std::size_t u = 0; u < kAlternativeTransformations.size(); ++u) {
    for (std::size_t o = 0; o < kAlternativeTransformations[u].size(); ++o) {
      auto const &src = kAlternativeTransformations[u][o];
      for (std::size_t i = 0; i < src.size(); ++i) {
        if (src[i] != 0) {
          out[u][o][i] = decode_packed(src[i]);
        }
      }
    }
  }
  return out;
}();

[[nodiscard]] MagneticSymmetryOperation
make_magnetic(DecodedMagneticOp const &d) noexcept {
  SymmetryOperation const op = make_operation(d.op);
  return {op.rotation, op.translation, d.time_reversal};
}

// Offset from a UNI number's smallest Hall setting to `hall_number`; the
// default (hall_number == 0) is offset 0. std::nullopt if the UNI number or the
// resulting offset is out of range.
[[nodiscard]] std::optional<int> hall_number_offset(int uni_number,
                                                    int hall_number) noexcept {
  if (uni_number < 1 || uni_number > kNumUniNumbers) {
    return std::nullopt;
  }
  const auto [num_halls, first_hall_number] =
      kMagneticUniMapping[static_cast<std::size_t>(uni_number)];
  int offset = 0;
  if (hall_number > 0 && hall_number <= kNumHallNumbers) {
    offset = hall_number - first_hall_number;
  } else if (hall_number != 0) {
    return std::nullopt;
  }
  if (offset < 0 || offset >= num_halls) {
    return std::nullopt;
  }
  return offset;
}

} // namespace

MagneticSymmetryOperations magnetic_operations_from_database(int uni_number,
                                                             int hall_number) {
  auto const offset = hall_number_offset(uni_number, hall_number);
  if (!offset) {
    return {};
  }
  auto const &idx =
      kMagneticOperationIndex[static_cast<std::size_t>(uni_number)]
                             [static_cast<std::size_t>(*offset)];
  auto const [count, start] = idx;
  auto const decoded =
      std::ranges::subrange(kDecodedMagneticOps.begin() + start,
                            kDecodedMagneticOps.begin() + start + count);
  MagneticSymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(count));
  std::ranges::transform(decoded, std::back_inserter(ops), make_magnetic);
  return ops;
}

SymmetryOperations magnetic_std_transformations(int uni_number,
                                                int hall_number) {
  auto const offset = hall_number_offset(uni_number, hall_number);
  if (!offset) {
    return {};
  }
  auto const &row =
      kDecodedAltTransformations[static_cast<std::size_t>(uni_number)]
                                [static_cast<std::size_t>(*offset)];

  // The identity transformation is always first; the rest run up to the
  // nullopt terminator.
  SymmetryOperations transforms;
  transforms.push_back({Matrix3i::Identity(), Vector3d::Zero()});
  for (auto const &slot : row | std::views::take_while([](auto const &s) {
                            return s.has_value();
                          })) {
    transforms.push_back(make_operation(*slot));
  }
  return transforms;
}

} // namespace cppcrystal::data
