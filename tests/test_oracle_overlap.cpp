#include "oracle.hpp"

#include <cppcrystal/core/overlap.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;

namespace {
Cell bcc(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0});
}

// Rutile-like tetragonal cell (lower symmetry, several operations).
Cell rutile() {
  Matrix3d lattice = Matrix3d::Zero();
  lattice(0, 0) = 4.6;
  lattice(1, 1) = 4.6;
  lattice(2, 2) = 3.0;
  Positions pos(6, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  pos.row(2) << 0.3, 0.3, 0.0;
  pos.row(3) << 0.7, 0.7, 0.0;
  pos.row(4) << 0.2, 0.8, 0.5;
  pos.row(5) << 0.8, 0.2, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0, 1, 1, 1, 1});
}
} // namespace

TEST_CASE("every reference symmetry op overlaps under our checker",
          "[oracle][overlap]") {
  for (Cell const &cell : {bcc(3.1), rutile()}) {
    OverlapChecker checker(cell, 1e-5);
    auto const ops = oracle::reference_symmetry(cell, 1e-5);
    REQUIRE(!ops.empty());
    for (auto const &op : ops) {
      INFO("rotation det = " << op.rotation.determinant());
      CHECK(checker.check_total_overlap(op.translation, op.rotation));
    }
  }
}

TEST_CASE("a corrupted reference op is rejected", "[oracle][overlap]") {
  Cell const cell = rutile();
  OverlapChecker checker(cell, 1e-5);
  auto ops = oracle::reference_symmetry(cell, 1e-5);
  REQUIRE(ops.size() > 1);
  // Shift a real operation's translation off by a clearly non-symmetric amount.
  SymmetryOperation bad = ops.front();
  bad.translation += Vector3d(0.13, 0.0, 0.0);
  CHECK_FALSE(checker.check_total_overlap(bad.translation, bad.rotation));
}
