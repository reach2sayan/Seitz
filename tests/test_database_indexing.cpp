// Equivalence test guarding the database refactor: the new map-based access
// (multi_index catalog + flat_map, via the public accessors) must agree with
// the original direct constexpr-array indexing for every Hall number. This is
// independent of the reference oracle, so it runs in the standalone suite.

#include <spglib/data/spacegroup_tables.hpp> // raw tables (old indexing)
#include <spglib/data/spg_database.hpp>      // accessors (new indexing)

#include <catch2/catch_test_macros.hpp>

#include <iterator>
#include <string_view>

using namespace spglib;

namespace {
// Old style: decode operations straight from the flat array using the
// {count, start} offsets.
SymmetryOperations operations_by_offset(int hall) {
  auto const &idx =
      data::kSymmetryOperationIndex[static_cast<std::size_t>(hall)];
  SymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(idx[0]));
  for (int i = 0; i < idx[0]; ++i)
    ops.push_back(data::decode_operation(
        data::kSymmetryOperations[static_cast<std::size_t>(idx[1] + i)]));
  return ops;
}
} // namespace

TEST_CASE("new map access matches old array indexing for all 530 Hall numbers",
          "[database]") {
  int op_mismatches = 0;
  int type_mismatches = 0;
  for (int hall = 1; hall <= data::kNumHallNumbers; ++hall) {
    // --- operations: offset-decoded vs flat_map ---
    auto const old_ops = operations_by_offset(hall);
    auto const &new_ops = data::database().operations.at(hall);
    bool ops_ok = old_ops.size() == new_ops.size();
    for (std::size_t s = 0; ops_ok && s < old_ops.size(); ++s)
      ops_ok = old_ops[s].rotation == new_ops[s].rotation &&
               old_ops[s].translation == new_ops[s].translation;
    op_mismatches += ops_ok ? 0 : 1;

    // --- metadata: raw array entry vs multi_index catalog ---
    auto const &raw = data::kSpacegroupTypes[static_cast<std::size_t>(hall)];
    auto const t = data::spacegroup_type(hall);
    bool const type_ok =
        t.hall_number == hall && t.number == raw.number &&
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
  auto const &by_number = data::database().catalog.get<data::ByNumber>();
  std::size_t total = 0;
  for (int n = 1; n <= 230; ++n) {
    auto const range = by_number.equal_range(n);
    total += static_cast<std::size_t>(std::distance(range.first, range.second));
  }
  CHECK(total == static_cast<std::size_t>(data::kNumHallNumbers));
  CHECK(std::distance(by_number.equal_range(3).first,
                      by_number.equal_range(3).second) == 3);
}
