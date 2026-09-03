// Self-consistency test for the Hall-symbol matcher. spglib's
// hal_match_hall_symbol_db is internal (no public oracle), so we exercise the
// matcher against the database itself: feeding a Hall setting's own database
// operations (the canonical setting) back through match_hall must
// recognise that Hall number, with a zero origin shift.

#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include "spacegroup/spacegroup.hpp"

#include "helpers.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <random>
#include <vector>

using namespace cppcrystal;
using cppcrystal::test::space_hall;

TEST_CASE(
    "each Hall setting matches its own database operations with zero shift",
    "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  int not_matched = 0;
  int nonzero_shift = 0;
  int first_unmatched = 0;
  for (int index = 1; index <= data::kNumHallNumbers; ++index) {
    HallNumber const hall = space_hall(index);
    auto const &ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    auto const shift =
        spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, centering, ops, 1e-5);
    if (!shift) {
      ++not_matched;
      if (first_unmatched == 0)
        first_unmatched = index;
      continue;
    }
    // Shift folded to [0,1) must be ~0 for the canonical setting.
    Vector3d const folded = math::wrap_to_unit_cell(*shift);
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
  for (int index : {1, 2, 3, 108, 349, 430, 462, 489, 523, 530}) {
    HallNumber const hall = space_hall(index);
    auto const &ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    INFO("hall " << index);
    CHECK(spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, centering, ops, 1e-5)
              .has_value());
  }
}

TEST_CASE("matching is independent of operation order", "[hall]") {
  // Exact database translations, shuffled: neither the rotation index nor the
  // origin-shift formula may care which operation carrying a generator
  // rotation comes first. The B-centered monoclinic / orthorhombic settings
  // (hall 13, 15, 34, ... 333) are the ones that bite: the database declares
  // them PRIMITIVE, so their B translation survives into dw and each
  // representative gives a different origin shift. No noise is added on
  // purpose — dw is folded to [0, 1) with a 1e-10 zero band (the reference's
  // mat_Dmod1), so a component that is 0 up to negative noise folds to ~1, a
  // representative the generator system need not admit.
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  std::mt19937 rng(99);
  std::vector<int> failures;
  for (int index = 1; index <= data::kNumHallNumbers; ++index) {
    HallNumber const hall = space_hall(index);
    Operations const db = data::operations_from_database(hall);
    std::vector<SymmetryOperation> shuffled(db.begin(), db.end());
    auto const centering = data::spacegroup_type(hall).centering;
    // Several orders per setting: one permutation can pick the good
    // representative by luck.
    for (int trial = 0; trial < 4; ++trial) {
      std::ranges::shuffle(shuffled, rng);
      if (!spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, centering, Operations{shuffled}, 1e-5)) {
        failures.push_back(index);
        break;
      }
    }
  }
  CAPTURE(failures);
  CHECK(failures.empty());
}

TEST_CASE("a corrupted operation is rejected", "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  for (int index : {3, 108, 349, 430, 462, 489, 523}) {
    HallNumber const hall = space_hall(index);
    Operations const db = data::operations_from_database(hall);
    std::vector<SymmetryOperation> ops(db.begin(), db.end());
    auto const centering = data::spacegroup_type(hall).centering;
    // Along all three axes: a component perpendicular to a rotation axis is
    // just a displaced axis (an origin shift, which the matcher rightly
    // accepts — for P2 that is any x or z), the parallel one is a screw.
    ops.back().translation += Vector3d(0.13, 0.13, 0.13);
    INFO("hall " << index);
    CHECK_FALSE(spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, centering, Operations{ops}, 1e-5));
  }
}

TEST_CASE("match_hall over the cubic settings", "[!benchmark]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  BENCHMARK("Ia-3d (hall 530, 96 ops)") {
    HallNumber const hall = space_hall(530);
    auto const &ops = data::operations_from_database(hall);
    return spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, data::spacegroup_type(hall).centering, ops, 1e-5);
  };
  BENCHMARK("Fm-3m (hall 523, 192 ops)") {
    HallNumber const hall = space_hall(523);
    auto const &ops = data::operations_from_database(hall);
    return spacegroup::SpacegroupMatcher<GroupFamily::space>::match_hall(
              bravais, hall, data::spacegroup_type(hall).centering, ops, 1e-5);
  };
}
