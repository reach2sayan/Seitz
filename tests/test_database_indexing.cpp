// Equivalence test guarding the database refactor: the constexpr catalog +
// lazily-built operations cache (via the public accessors) must agree with the
// original direct constexpr-array indexing for every Hall number. This is
// independent of the reference oracle, so it runs in the standalone suite.

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/data/spacegroup_metadata_tables.hpp>  // raw metadata (old indexing)
#include "data/spacegroup_operation_tables.hpp" // raw ops (old indexing)
#include <cppcrystal/data/spg_database.hpp>                // accessors (new indexing)

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

#include <string_view>

using namespace cppcrystal;
using cppcrystal::test::space_hall;

namespace {
// Reference unpack of one packed operation (spgdb_decode_symmetry: rotation
// base-3 {-1,0,1}, translation base-12 n/12). Kept local to this test so the
// equivalence check has an independent decoder to compare the database against.
SymmetryOperation decode_operation(int encoded) {
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

// Old style: decode operations straight from the flat array using the
// {count, start} offsets.
Operations operations_by_offset(int hall) {
  auto const &idx =
      data::kSymmetryOperationIndex[static_cast<std::size_t>(hall)];
  std::vector<SymmetryOperation> ops;
  ops.reserve(static_cast<std::size_t>(idx[0]));
  for (int i = 0; i < idx[0]; ++i)
    ops.push_back(decode_operation(
        data::kSymmetryOperations[static_cast<std::size_t>(idx[1] + i)]));
  return Operations{std::move(ops)};
}
} // namespace

TEST_CASE("new map access matches old array indexing for all 530 Hall numbers",
          "[database]") {
  int op_mismatches = 0;
  int type_mismatches = 0;
  for (int index = 1; index <= data::kNumHallNumbers; ++index) {
    HallNumber const hall = space_hall(index);
    // --- operations: offset-decoded vs lazily-built cache ---
    auto const old_ops = operations_by_offset(index);
    auto const &new_ops = data::operations_from_database(hall);
    bool ops_ok = old_ops.size() == new_ops.size();
    for (std::size_t s = 0; ops_ok && s < old_ops.size(); ++s)
      ops_ok = old_ops[s].rotation == new_ops[s].rotation &&
               old_ops[s].translation == new_ops[s].translation;
    op_mismatches += ops_ok ? 0 : 1;

    // --- metadata: raw array entry vs constexpr catalog ---
    auto const &raw = data::kSpacegroupTypes[static_cast<std::size_t>(index)];
    auto const &t = data::spacegroup_type(hall);
    bool const type_ok =
        t.number == raw.number &&
        t.schoenflies == std::string_view(raw.schoenflies) &&
        t.hall_symbol == std::string_view(raw.hall_symbol) &&
        t.international == std::string_view(raw.international) &&
        t.international_full == std::string_view(raw.international_full) &&
        t.international_short == std::string_view(raw.international_short) &&
        t.choice == std::string_view(raw.choice) &&
        static_cast<int>(t.centering) == raw.centering &&
        t.pointgroup_number == raw.pointgroup_number;
    type_mismatches += type_ok ? 0 : 1;
  }
  CHECK(op_mismatches == 0);
  CHECK(type_mismatches == 0);
}

TEST_CASE("the by-number catalog index enumerates every Hall setting",
          "[database]") {
  // Total settings across all 230 space-group numbers must be exactly 530, and
  // space group 3 (P2) has three Hall settings (b, c, a unique-axis choices).
  std::size_t total = 0;
  for (int n = 1; n <= 230; ++n) {
    total += data::halls_with_number<GroupFamily::space>(n).size();
  }
  CHECK(total == static_cast<std::size_t>(data::kNumHallNumbers));
  CHECK(data::halls_with_number<GroupFamily::space>(3).size() == 3);

  // The catalog is a genuine constant expression.
  static_assert(data::spacegroup_type(space_hall(1)).number == 1);
  static_assert(data::halls_with_number<GroupFamily::space>(3).size() == 3);
}
