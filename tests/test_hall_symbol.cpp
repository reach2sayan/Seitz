// Self-consistency test for the Hall-symbol matcher. spglib's
// hal_match_hall_symbol_db is internal (no public oracle), so we exercise the
// matcher against the database itself: feeding a Hall setting's own database
// operations (the canonical setting) back through match_hall_symbol must
// recognise that Hall number, with a zero origin shift.

#include <spglib/data/spg_database.hpp>
#include <spglib/spacegroup/hall_symbol.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace spglib;

TEST_CASE(
    "each Hall setting matches its own database operations with zero shift",
    "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  int not_matched = 0;
  int nonzero_shift = 0;
  int first_unmatched = 0;
  for (int hall = 1; hall <= data::kNumHallNumbers; ++hall) {
    auto const &ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    auto const shift =
        spacegroup::match_hall_symbol(bravais, hall, centering, ops, 1e-5);
    if (!shift) {
      ++not_matched;
      if (first_unmatched == 0)
        first_unmatched = hall;
      continue;
    }
    // Shift folded to [0,1) must be ~0 for the canonical setting.
    Vector3d const folded = math::mod1(*shift);
    if (folded.cwiseAbs().maxCoeff() > 1e-6 &&
        (folded.array() - 1.0).abs().maxCoeff() > 1e-6)
      ++nonzero_shift;
  }
  INFO("first unmatched hall = "
       << first_unmatched << ", nonzero-shift count = " << nonzero_shift);
  CHECK(not_matched == 0);
  CHECK(nonzero_shift == 0);
}

TEST_CASE("a representative selection of Hall settings match", "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 4.0;
  for (int hall : {1, 2, 3, 108, 349, 430, 462, 489, 523, 530}) {
    auto const &ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    INFO("hall " << hall);
    CHECK(spacegroup::match_hall_symbol(bravais, hall, centering, ops, 1e-5)
              .has_value());
  }
}
