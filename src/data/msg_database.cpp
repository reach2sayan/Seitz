#include <cppcrystal/data/msg_database.hpp>

#include <cppcrystal/core/mdspan.hpp>

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
    out[k] = {.op = decode_packed(enc % kTimeReversalBase),
              .time_reversal = enc / kTimeReversalBase != 0};
  }
  return out;
}();

// Every alternative standardized-setting transformation, decoded once at
// compile time into a flat [uni][offset][slot] table. std::nullopt is the
// per-row terminator (the source table zero-terminates each [uni][offset]
// row); using optional rather than a magic 0 keeps the no-sentinel house rule.
constexpr std::size_t kNumAltSlots = 7;
constexpr std::size_t kMaxHallsPerUni = 18;
using AltTable = md::table<std::optional<DecodedOp>, kNumUniNumbers + 1,
                           kMaxHallsPerUni, kNumAltSlots>;
constexpr auto kDecodedAltTransformations = [] {
  std::array<std::optional<DecodedOp>,
             (kNumUniNumbers + 1) * kMaxHallsPerUni * kNumAltSlots>
      out{};
  md::mdspan<std::optional<DecodedOp>, AltTable::extents_type> const view(
      out.data());
  for (std::size_t u = 0; u < kAlternativeTransformations.size(); ++u) {
    for (std::size_t o = 0; o < kAlternativeTransformations[u].size(); ++o) {
      auto const &src = kAlternativeTransformations[u][o];
      for (std::size_t i = 0; i < src.size(); ++i) {
        if (src[i] != 0) {
          view[u, o, i] = decode_packed(src[i]);
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
// default (hall_number == 0) is offset 0. std::nullopt if the UNI number, the
// Hall number or the resulting offset is out of range.
[[nodiscard]] std::optional<int> hall_number_offset(int uni_number,
                                                    int hall_number) noexcept {
  if (uni_number < 1 || uni_number > kNumUniNumbers) {
    return std::nullopt;
  }
  if (hall_number != 0 && (hall_number < 1 || hall_number > kNumHallNumbers)) {
    return std::nullopt;
  }
  const auto [num_halls, first_hall_number] =
      kMagneticUniMapping[static_cast<std::size_t>(uni_number)];
  int const offset = hall_number == 0 ? 0 : hall_number - first_hall_number;
  if (offset < 0 || offset >= num_halls) {
    return std::nullopt;
  }
  return offset;
}

// Flat (uni, offset) slot of a per-setting cache; slot 0 (uni 0) is never
// filled and serves as the empty fallback.
[[nodiscard]] std::size_t setting_slot(int uni_number, int offset) noexcept {
  return static_cast<std::size_t>(uni_number) * kMaxHallsPerUni +
         static_cast<std::size_t>(offset);
}

// Materialise `build(uni, offset)` for every valid setting once.
template <class T, class Build>
[[nodiscard]] std::vector<T> build_setting_table(Build build) {
  std::vector<T> table((kNumUniNumbers + 1) * kMaxHallsPerUni);
  for (int uni = 1; uni <= kNumUniNumbers; ++uni) {
    int const num_halls = kMagneticUniMapping[static_cast<std::size_t>(uni)][0];
    for (int offset = 0; offset < num_halls; ++offset) {
      table[setting_slot(uni, offset)] = build(uni, offset);
    }
  }
  return table;
}

[[nodiscard]] MagneticSymmetryOperations build_operations(int uni_number,
                                                          int offset) {
  auto const &idx = kMagneticOperationIndex[static_cast<std::size_t>(uni_number)]
                                           [static_cast<std::size_t>(offset)];
  auto const [count, start] = idx;
  auto const decoded =
      std::ranges::subrange(kDecodedMagneticOps.begin() + start,
                            kDecodedMagneticOps.begin() + start + count);
  MagneticSymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(count));
  std::ranges::transform(decoded, std::back_inserter(ops), make_magnetic);
  return ops;
}

[[nodiscard]] SymmetryOperations build_transformations(int uni_number,
                                                       int offset) {
  AltTable const table(kDecodedAltTransformations.data());
  auto const row = md::submdspan(table, static_cast<std::size_t>(uni_number),
                                 static_cast<std::size_t>(offset),
                                 md::full_extent);
  // The identity transformation is always first; the rest run up to the
  // nullopt terminator.
  SymmetryOperations transforms;
  transforms.push_back(
      {.rotation = Matrix3i::Identity(), .translation = Vector3d::Zero()});
  for (std::size_t i = 0; i < row.extent(0) && row[i]; ++i) {
    transforms.push_back(make_operation(*row[i]));
  }
  return transforms;
}

} // namespace

MagneticSymmetryOperations const &
magnetic_operations_from_database(int uni_number, int hall_number) {
  static auto const table =
      build_setting_table<MagneticSymmetryOperations>(build_operations);
  auto const offset = hall_number_offset(uni_number, hall_number);
  return table[offset ? setting_slot(uni_number, *offset) : 0];
}

SymmetryOperations const &magnetic_std_transformations(int uni_number,
                                                       int hall_number) {
  static auto const table =
      build_setting_table<SymmetryOperations>(build_transformations);
  auto const offset = hall_number_offset(uni_number, hall_number);
  return table[offset ? setting_slot(uni_number, *offset) : 0];
}

} // namespace cppcrystal::data
