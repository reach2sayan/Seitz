#include "oracle.hpp"

#include <spglib/symmetry/find_symmetry.hpp>
#include <spglib/symmetry/pointgroup.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace spglib;

namespace {
Cell primitive_cubic(double a) {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Matrix3d::Identity() * a, pos, {0});
}

Cell rutile() {
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.6;
  l(1, 1) = 4.6;
  l(2, 2) = 3.0;
  Positions pos(6, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  pos.row(2) << 0.3, 0.3, 0.0;
  pos.row(3) << 0.7, 0.7, 0.0;
  pos.row(4) << 0.8, 0.2, 0.5;
  pos.row(5) << 0.2, 0.8, 0.5;
  return Cell(l, pos, {0, 0, 1, 1, 1, 1});
}

Cell hexagonal() {
  Matrix3d l;
  l.col(0) = Vector3d(3.0, 0.0, 0.0);
  l.col(1) = Vector3d(-1.5, 2.598076211, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, 5.0);
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(l, pos, {0});
}
} // namespace

TEST_CASE("get_pointgroup matches spg_get_pointgroup", "[oracle][pointgroup]") {
  for (Cell const &cell : {primitive_cubic(4.0), rutile(), hexagonal()}) {
    auto ops = symmetry::find_symmetry(cell, 1e-5);
    REQUIRE(ops);
    // Feed identical rotations to both implementations so axis selection and
    // de-duplication order match.
    auto const rotations = symmetry::rotations_of(*ops);
    auto ours = symmetry::get_pointgroup(rotations);
    REQUIRE(ours);
    auto const ref = oracle::reference_pointgroup(rotations);

    INFO("ours = " << ours->pointgroup.symbol << " (" << ours->pointgroup.number
                   << "), ref = " << ref.symbol << " (" << ref.number << ")");
    CHECK(ours->pointgroup.number == ref.number);
    CHECK(std::string(ours->pointgroup.symbol) == ref.symbol);
    CHECK(ours->transformation == ref.transformation);
  }
}
