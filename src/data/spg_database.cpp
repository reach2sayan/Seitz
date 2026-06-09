#include <spglib/data/spg_database.hpp>

#include <spglib/data/spacegroup_tables.hpp>

#include <cassert>
#include <utility>

namespace spglib::data {

SymmetryOperation decode_operation(int encoded) noexcept {
  Matrix3i rot;
  int r = encoded % 19683; // 3^9
  int digit = 6561;        // 3^8
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      rot(i, j) = (r % (digit * 3)) / digit - 1;
      digit /= 3;
    }

  int const t = encoded / 19683;
  Vector3d trans;
  digit = 144; // 12^2
  for (int i = 0; i < 3; ++i) {
    trans[i] = static_cast<double>((t % (digit * 12)) / digit) / 12.0;
    digit /= 12;
  }
  return {rot, trans};
}

namespace {

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
    for (int i = 0; i < idx[0]; ++i)
      ops.push_back(decode_operation(
          kSymmetryOperations[static_cast<std::size_t>(idx[1] + i)]));
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
