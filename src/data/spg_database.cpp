#include <spglib/data/spg_database.hpp>

#include <spglib/data/spacegroup_tables.hpp>

#include <array>
#include <cassert>
#include <cstdint>
#include <ranges>
#include <utility>

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
    int r = encoded % 19683; // 3^9
    int digit = 6561;        // 3^8
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

  const Vector3d trans{d.trans_num[0] / 12.0, d.trans_num[1] / 12.0,
                       d.trans_num[2] / 12.0};

  return {rot, trans};
}

[[nodiscard]] Database build_database() {
  Database db;
  for (int hall = 1; hall <= kNumHallNumbers; ++hall) {
    auto const &r = kSpacegroupTypes[static_cast<std::size_t>(hall)];
    db.catalog.insert(SpacegroupType{
        hall, r.number, r.schoenflies, r.hall_symbol, r.international,
        r.international_full, r.international_short, r.choice,
        static_cast<Centering>(r.centering), r.pointgroup_number});

    auto const &idx = kSymmetryOperationIndex[static_cast<std::size_t>(hall)];
    SymmetryOperations ops;
    ops.reserve(static_cast<std::size_t>(idx[0]));

    std::ranges::transform(
        std::views::iota(0, idx[0]), std::back_inserter(ops), [&](int i) {
          return make_operation(
              kDecodedOps[static_cast<std::size_t>(idx[1] + i)]);
        });
    db.operations.emplace(hall, std::move(ops));
  }
  // Integrity check (debug only): every Hall setting must be present once. A
  // duplicate key would have made a hashed_unique / flat_map insert silently
  // fail, shrinking these below 530.
  assert(db.catalog.size() == static_cast<std::size_t>(kNumHallNumbers));
  assert(db.operations.size() == static_cast<std::size_t>(kNumHallNumbers));
  return db;
}

} // namespace

Database const &database() {
  static Database const db = build_database();
  return db;
}

SymmetryOperations operations_from_database(int hall_number) {
  auto const &ops = database().operations;
  auto const it = ops.find(hall_number);
  return it == ops.end() ? SymmetryOperations{} : it->second;
}

SpacegroupType spacegroup_type(int hall_number) {
  auto const &by_hall = database().catalog.get<ByHall>();
  auto const it = by_hall.find(hall_number);
  return it == by_hall.end() ? SpacegroupType{} : *it;
}

} // namespace spglib::data
