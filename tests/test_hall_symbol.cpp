// Self-consistency test for the Hall-symbol matcher. spglib's
// hal_match_hall_symbol_db is internal (no public oracle), so we exercise the
// matcher against the database itself: feeding a Hall setting's own database
// operations (the canonical setting) back through match_hall_symbol must
// recognise that Hall number, with a zero origin shift.

#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/spacegroup/hall_symbol.hpp>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <random>

using namespace cppcrystal;

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
  for (int hall : {1, 2, 3, 108, 349, 430, 462, 489, 523, 530}) {
    auto const &ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    INFO("hall " << hall);
    CHECK(spacegroup::match_hall_symbol(bravais, hall, centering, ops, 1e-5)
              .has_value());
  }
}

TEST_CASE("matching is independent of operation order and survives "
          "sub-tolerance noise",
          "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  std::mt19937 rng(99);
  // symprec is Cartesian (1e-5 on a 3 A cell = 3.3e-6 fractional) and the
  // origin shift VSpU . dw amplifies generator noise by a small factor, so the
  // per-component noise has to sit well below that.
  std::uniform_real_distribution<double> noise(-1e-7, 1e-7);
  int failures = 0;
  for (int hall = 1; hall <= data::kNumHallNumbers; ++hall) {
    SymmetryOperations ops = data::operations_from_database(hall);
    std::ranges::shuffle(ops, rng);
    for (auto &op : ops) {
      op.translation += Vector3d(noise(rng), noise(rng), noise(rng));
    }
    auto const centering = data::spacegroup_type(hall).centering;
    if (!spacegroup::match_hall_symbol(bravais, hall, centering, ops, 1e-5)) {
      ++failures;
    }
  }
  CHECK(failures == 0);
}

TEST_CASE("a corrupted operation is rejected", "[hall]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  for (int hall : {3, 108, 349, 430, 462, 489, 523}) {
    SymmetryOperations ops = data::operations_from_database(hall);
    auto const centering = data::spacegroup_type(hall).centering;
    // Along all three axes: a component perpendicular to a rotation axis is
    // just a displaced axis (an origin shift, which the matcher rightly
    // accepts — for P2 that is any x or z), the parallel one is a screw.
    ops.back().translation += Vector3d(0.13, 0.13, 0.13);
    INFO("hall " << hall);
    CHECK_FALSE(
        spacegroup::match_hall_symbol(bravais, hall, centering, ops, 1e-5));
  }
}

TEST_CASE("match_hall_symbol over the cubic settings", "[!benchmark]") {
  Matrix3d const bravais = Matrix3d::Identity() * 3.0;
  BENCHMARK("Ia-3d (hall 530, 96 ops)") {
    auto const &ops = data::operations_from_database(530);
    return spacegroup::match_hall_symbol(
        bravais, 530, data::spacegroup_type(530).centering, ops, 1e-5);
  };
  BENCHMARK("Fm-3m (hall 523, 192 ops)") {
    auto const &ops = data::operations_from_database(523);
    return spacegroup::match_hall_symbol(
        bravais, 523, data::spacegroup_type(523).centering, ops, 1e-5);
  };
}
