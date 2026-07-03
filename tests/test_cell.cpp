#include <cppcrystal/core/cell.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using Catch::Approx;

namespace {
// Simple cubic cell, lattice parameter a, with one atom at the origin.
Cell cubic(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(lattice, pos, {0});
}
} // namespace

TEST_CASE("Cell exposes size, types, positions", "[cell]") {
  Matrix3d lattice = Matrix3d::Identity() * 3.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  Cell c(lattice, pos, {0, 1});
  CHECK(c.size() == 2);
  CHECK(c.type(1) == 1);
  CHECK(c.position(1).isApprox(Vector3d(0.5, 0.5, 0.5)));
  CHECK_FALSE(c.aperiodic_axis().has_value());
}

TEST_CASE("Cell volume is |det(lattice)|", "[cell]") {
  Cell c = cubic(2.0);
  CHECK(c.volume() == Approx(8.0));
  Matrix3d l;
  l << 0, 1, 0, -1, 0, 0, 0, 0, 1; // det = 1, still volume 1
  Cell c2(l, Positions(0, 3), {});
  CHECK(c2.volume() == Approx(1.0));
}

TEST_CASE("Cell to_cartesian multiplies by the lattice", "[cell]") {
  Matrix3d lattice = Matrix3d::Identity() * 4.0;
  Cell c(lattice, Positions(0, 3), {});
  CHECK(c.to_cartesian(Vector3d(0.5, 0.25, 0.0))
            .isApprox(Vector3d(2.0, 1.0, 0.0)));
}

TEST_CASE("Cell metric tensor is diagonal for an orthogonal lattice",
          "[cell]") {
  Matrix3d lattice;
  lattice << 2, 0, 0, 0, 3, 0, 0, 0, 4;
  Matrix3d g = Cell(lattice, Positions(0, 3), {}).metric();
  CHECK(g(0, 0) == Approx(4.0));
  CHECK(g(1, 1) == Approx(9.0));
  CHECK(g(2, 2) == Approx(16.0));
  CHECK(g(0, 1) == Approx(0.0));
}
