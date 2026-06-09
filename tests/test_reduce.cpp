#include <spglib/reduce/delaunay.hpp>
#include <spglib/reduce/niggli.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace spglib;
using Catch::Approx;

namespace {
// A sheared (non-reduced) triclinic lattice, columns = basis vectors.
Matrix3d sheared() {
  Matrix3d l;
  l.col(0) = Vector3d(1.0, 0.0, 0.0);
  l.col(1) = Vector3d(2.0, 1.0, 0.0); // strongly sheared
  l.col(2) = Vector3d(3.0, 1.0, 1.0);
  return l;
}
} // namespace

TEST_CASE("delaunay reduction preserves cell volume", "[reduce]") {
  auto red = reduce::delaunay_reduce(sheared(), 1e-5);
  REQUIRE(red);
  CHECK(std::abs(red->determinant()) ==
        Approx(std::abs(sheared().determinant())));
}

TEST_CASE("delaunay reduction yields a right-handed cell", "[reduce]") {
  auto red = reduce::delaunay_reduce(sheared(), 1e-5);
  REQUIRE(red);
  CHECK(red->determinant() > 0.0);
}

TEST_CASE("delaunay reduction shortens the basis", "[reduce]") {
  Matrix3d const in = sheared();
  auto red = reduce::delaunay_reduce(in, 1e-5);
  REQUIRE(red);
  double const in_len = in.colwise().norm().sum();
  double const out_len = red->colwise().norm().sum();
  CHECK(out_len <= in_len + 1e-9);
}

TEST_CASE("niggli reduction preserves cell volume", "[reduce]") {
  auto red = reduce::niggli_reduce(sheared(), 1e-5);
  REQUIRE(red);
  CHECK(std::abs(red->determinant()) ==
        Approx(std::abs(sheared().determinant())));
}

TEST_CASE("niggli reduced cell satisfies A <= B <= C", "[reduce]") {
  auto red = reduce::niggli_reduce(sheared(), 1e-5);
  REQUIRE(red);
  Matrix3d const g = red->transpose() * (*red);
  double const eps = 1e-8;
  CHECK(g(0, 0) <= g(1, 1) + eps);
  CHECK(g(1, 1) <= g(2, 2) + eps);
}

TEST_CASE("reducing an already-cubic lattice is a no-op in length",
          "[reduce]") {
  Matrix3d const cube = Matrix3d::Identity() * 4.0;
  auto d = reduce::delaunay_reduce(cube, 1e-5);
  auto n = reduce::niggli_reduce(cube, 1e-5);
  REQUIRE(d);
  REQUIRE(n);
  CHECK(std::abs(d->determinant()) == Approx(64.0));
  CHECK(std::abs(n->determinant()) == Approx(64.0));
}
