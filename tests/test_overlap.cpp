#include <cppcrystal/core/overlap.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;

namespace {
// Body-centered cubic: atoms at (0,0,0) and (1/2,1/2,1/2).
Cell bcc(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(lattice, pos, {0, 0});
}

Matrix3i inversion() { return -Matrix3i::Identity(); }
} // namespace

TEST_CASE("is_overlap uses the Cartesian minimal image", "[overlap]") {
  Matrix3d lat = Matrix3d::Identity() * 5.0;
  CHECK(
      is_overlap(Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, -2.0, 3.0), lat, 1e-5));
  CHECK_FALSE(
      is_overlap(Vector3d(0.0, 0.0, 0.0), Vector3d(0.1, 0.0, 0.0), lat, 1e-5));
  // 0.1 fractional along a=5 -> 0.5 Angstrom > symprec.
}

TEST_CASE("identity is always a symmetry", "[overlap]") {
  OverlapChecker checker(bcc(3.0));
  CHECK(checker.check_total_overlap(Vector3d::Zero(), Matrix3i::Identity(),
                                    1e-5));
}

TEST_CASE("inversion is a symmetry of bcc", "[overlap]") {
  OverlapChecker checker(bcc(3.0));
  CHECK(checker.check_total_overlap(Vector3d::Zero(), inversion(), 1e-5));
}

TEST_CASE("a bogus translation is rejected", "[overlap]") {
  OverlapChecker checker(bcc(3.0));
  CHECK_FALSE(checker.check_total_overlap(Vector3d(0.3, 0.0, 0.0),
                                          Matrix3i::Identity(), 1e-5));
}

TEST_CASE("the bcc centering translation is a symmetry", "[overlap]") {
  // (1/2,1/2,1/2) maps atom 0 -> atom 1 and atom 1 -> atom 0 (mod lattice).
  OverlapChecker checker(bcc(3.0));
  CHECK(checker.check_total_overlap(Vector3d(0.5, 0.5, 0.5),
                                    Matrix3i::Identity(), 1e-5));
}

TEST_CASE("symmetry breaks when the two atoms differ in type", "[overlap]") {
  Matrix3d lattice = Matrix3d::Identity() * 3.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  OverlapChecker checker(
      Cell(lattice, pos, {0, 1})); // CsCl-like, distinct types
  // The centering translation now swaps unlike atoms -> not a symmetry.
  CHECK_FALSE(checker.check_total_overlap(Vector3d(0.5, 0.5, 0.5),
                                          Matrix3i::Identity(), 1e-5));
  // Inversion keeps each atom in place -> still a symmetry.
  CHECK(checker.check_total_overlap(Vector3d::Zero(), -Matrix3i::Identity(),
                                    1e-5));
}
