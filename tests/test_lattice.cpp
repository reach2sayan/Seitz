#include <cppcrystal/core/lattice.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using Catch::Approx;

namespace {
// A sheared (non-reduced) triclinic lattice, columns = basis vectors.
Lattice sheared() {
  Matrix3d l;
  l.col(0) = Vector3d(1.0, 0.0, 0.0);
  l.col(1) = Vector3d(2.0, 1.0, 0.0); // strongly sheared
  l.col(2) = Vector3d(3.0, 1.0, 1.0);
  return Lattice{l};
}
} // namespace

TEST_CASE("from_basis rejects a singular basis", "[lattice]") {
  CHECK_FALSE(Lattice::from_basis(Matrix3d::Zero()).has_value());
  Matrix3d degenerate = Matrix3d::Identity();
  degenerate.col(2) = degenerate.col(1); // two identical basis vectors
  CHECK_FALSE(Lattice::from_basis(degenerate).has_value());
  REQUIRE(Lattice::from_basis(Matrix3d::Identity() * 3.0).has_value());
}

TEST_CASE("volume is |det(basis)|", "[lattice]") {
  CHECK(Lattice{Matrix3d::Identity() * 2.0}.volume() == Approx(8.0));
  Matrix3d l;
  l << 0, 1, 0, -1, 0, 0, 0, 0, 1; // det = 1, still volume 1
  CHECK(Lattice{l}.volume() == Approx(1.0));
}

TEST_CASE("metric tensor is L^T . L", "[lattice]") {
  Matrix3d l;
  l << 2, 0, 0, 0, 3, 0, 0, 0, 4;
  Matrix3d const g = Lattice{l}.metric();
  CHECK(g(0, 0) == Approx(4.0));
  CHECK(g(1, 1) == Approx(9.0));
  CHECK(g(2, 2) == Approx(16.0));
  CHECK(g(0, 1) == Approx(0.0));
}

TEST_CASE("cartesian and fractional coordinates round-trip", "[lattice]") {
  Lattice const l{Matrix3d::Identity() * 4.0};
  CHECK(l.to_cartesian(Vector3d(0.5, 0.25, 0.0))
            .isApprox(Vector3d(2.0, 1.0, 0.0)));
  Vector3d const frac(0.3, -0.7, 1.4);
  CHECK(l.to_fractional(l.to_cartesian(frac)).isApprox(frac, 1e-12));
}

TEST_CASE("transformed right-multiplies the basis", "[lattice]") {
  Lattice const l = sheared();
  Matrix3d t;
  t << 0, 1, 0, 1, 0, 0, 0, 0, 1; // swap the first two basis vectors
  CHECK(l.transformed(t).matrix().col(0).isApprox(l.matrix().col(1)));
  CHECK(l.transformed(Matrix3d::Identity()).matrix().isApprox(l.matrix()));
}

TEST_CASE("rigid_rotation_to maps one lattice frame onto the other",
          "[lattice]") {
  Lattice const from{Matrix3d::Identity() * 2.0};
  // 90 degrees about z.
  Matrix3d rot;
  rot << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  Lattice const to{rot * from.matrix()};
  CHECK(from.rigid_rotation_to(to).isApprox(rot, 1e-12));
  CHECK(from.rigid_rotation_to(from).isApprox(Matrix3d::Identity(), 1e-12));
}

TEST_CASE("delaunay reduction preserves cell volume", "[lattice][reduce]") {
  auto red = sheared().delaunay(1e-5);
  REQUIRE(red);
  CHECK(red->volume() == Approx(sheared().volume()));
}

TEST_CASE("delaunay reduction yields a right-handed cell",
          "[lattice][reduce]") {
  auto red = sheared().delaunay(1e-5);
  REQUIRE(red);
  CHECK(red->matrix().determinant() > 0.0);
}

TEST_CASE("delaunay reduction shortens the basis", "[lattice][reduce]") {
  Lattice const in = sheared();
  auto red = in.delaunay(1e-5);
  REQUIRE(red);
  double const in_len = in.matrix().colwise().norm().sum();
  double const out_len = red->matrix().colwise().norm().sum();
  CHECK(out_len <= in_len + 1e-9);
}

TEST_CASE("2D delaunay reduction leaves the unique axis spanning the same line",
          "[lattice][reduce]") {
  Lattice const in = sheared();
  auto red = in.delaunay_in_plane(2, 1e-5);
  REQUIRE(red);
  // The unique axis is kept (up to the sign flip that keeps the cell
  // right-handed); the other two are reduced within their plane.
  CHECK(
      red->matrix().col(2).cwiseAbs().isApprox(in.matrix().col(2).cwiseAbs()));
  CHECK(red->volume() == Approx(in.volume()));
}

TEST_CASE("niggli reduction preserves cell volume", "[lattice][reduce]") {
  auto red = sheared().niggli(1e-5);
  REQUIRE(red);
  CHECK(red->volume() == Approx(sheared().volume()));
}

TEST_CASE("niggli reduced cell satisfies A <= B <= C", "[lattice][reduce]") {
  auto red = sheared().niggli(1e-5);
  REQUIRE(red);
  Matrix3d const g = red->metric();
  double const eps = 1e-8;
  CHECK(g(0, 0) <= g(1, 1) + eps);
  CHECK(g(1, 1) <= g(2, 2) + eps);
}

TEST_CASE("reducing an already-cubic lattice is a no-op in length",
          "[lattice][reduce]") {
  Lattice const cube{Matrix3d::Identity() * 4.0};
  auto d = cube.delaunay(1e-5);
  auto n = cube.niggli(1e-5);
  REQUIRE(d);
  REQUIRE(n);
  CHECK(d->volume() == Approx(64.0));
  CHECK(n->volume() == Approx(64.0));
}
