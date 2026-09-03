#include "oracle.hpp"
#include "helpers.hpp"

#include <cppcrystal/data/spg_database.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using cppcrystal::test::space_hall;

TEST_CASE("db ops match spg_get_symmetry_from_database for all 530 Hall numbers",
          "[oracle][database]") {
  int mismatches = 0;
  int first_bad = 0;
  for (int index = 1; index <= data::kNumHallNumbers; ++index) {
    auto const &ours = data::operations_from_database(space_hall(index));
    auto const ref = oracle::reference_database_operations(index);
    bool ok = ours.size() == ref.size();
    for (std::size_t s = 0; ok && s < ours.size(); ++s)
      ok = ours[s].rotation == ref[s].rotation &&
           (ours[s].translation - ref[s].translation).cwiseAbs().maxCoeff() <
               1e-12;
    if (!ok && first_bad == 0)
      first_bad = index;
    mismatches += ok ? 0 : 1;
  }
  INFO("first mismatching hall number = " << first_bad);
  CHECK(mismatches == 0);
}

TEST_CASE("db spg types match spg_get_spacegroup_type for all 530 Hall numbers",
          "[oracle][database]") {
  int mismatches = 0;
  int first_bad = 0;
  for (int index = 1; index <= data::kNumHallNumbers; ++index) {
    auto const &ours = data::spacegroup_type(space_hall(index));
    auto const ref = oracle::reference_spacegroup_type(index);
    bool const ok =
        ours.number == ref.number && ours.hall_symbol == ref.hall_symbol &&
        ours.international == ref.international &&
        ours.international_full == ref.international_full &&
        ours.international_short == ref.international_short &&
        ours.schoenflies == ref.schoenflies && ours.choice == ref.choice;
    if (!ok && first_bad == 0)
      first_bad = index;
    mismatches += ok ? 0 : 1;
  }
  INFO("first mismatching hall number = " << first_bad);
  CHECK(mismatches == 0);
}
