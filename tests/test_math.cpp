#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using Catch::Approx;

TEST_CASE("nint rounds half away from zero", "[math]") {
  CHECK(math::nint(2.4) == 2);
  CHECK(math::nint(2.5) == 3);
  CHECK(math::nint(2.6) == 3);
  CHECK(math::nint(-2.4) == -2);
  CHECK(math::nint(-2.5) == -3);
  CHECK(math::nint(0.5) == 1);
  CHECK(math::nint(-0.5) == -1);
  CHECK(math::nint(0.0) == 0);
}

TEST_CASE("wrap_to_unit_cell folds a coordinate into the unit interval",
          "[math]") {
  CHECK(math::wrap_to_unit_cell(0.2) == Approx(0.2));
  CHECK(math::wrap_to_unit_cell(1.2) == Approx(0.2));
  CHECK(math::wrap_to_unit_cell(-0.2) == Approx(0.8));
  CHECK(math::wrap_to_unit_cell(2.0) == Approx(0.0));
  CHECK(math::wrap_to_unit_cell(-1.0) == Approx(0.0));
  // Just below zero stays near zero, not near one.
  CHECK(math::wrap_to_unit_cell(-1e-11) == Approx(0.0).margin(1e-9));
}

TEST_CASE("nearest_offset returns the centered remainder x - nint(x)",
          "[math]") {
  Vector3d const off = math::nearest_offset(Vector3d(0.3, 0.7, 3.25));
  CHECK(off[0] == Approx(0.3));
  CHECK(off[1] == Approx(-0.3));
  CHECK(off[2] == Approx(0.25));
  CHECK(math::nearest_offset(Vector3d(-0.7, 0.0, 0.0))[0] == Approx(0.3));
}

TEST_CASE("is_int_matrix respects symprec", "[math]") {
  CHECK(math::is_int_matrix(Matrix3d::Identity(), 1e-5));
  Matrix3d m = Matrix3d::Identity();
  m(0, 1) = 0.5;
  CHECK_FALSE(math::is_int_matrix(m, 1e-5));
  m(0, 1) = 1.0 + 1e-7;
  CHECK(math::is_int_matrix(m, 1e-5));
}

TEST_CASE("round_to_int casts element-wise via nint", "[math]") {
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = 0.9;
  m(1, 1) = 2.4;
  m(2, 2) = -1.5;
  Matrix3i r = math::round_to_int(m);
  CHECK(r(0, 0) == 1);
  CHECK(r(1, 1) == 2);
  CHECK(r(2, 2) == -2);
}

TEST_CASE("checked inverse fails on a singular matrix", "[math]") {
  CHECK_FALSE(math::inverse(Matrix3d::Zero(), 1e-10).has_value());
  Matrix3d a;
  a << 2, 0, 0, 0, 4, 0, 0, 0, 8;
  auto inv = math::inverse(a, 1e-10);
  REQUIRE(inv.has_value());
  CHECK((*inv * a).isApprox(Matrix3d::Identity(), 1e-12));
}

TEST_CASE("metric tensor is L^T . L", "[math]") {
  Matrix3d lattice;
  lattice << 2, 0, 0, 0, 3, 0, 0, 0, 4; // columns = basis vectors
  Matrix3d g = math::metric_tensor(lattice);
  CHECK(g(0, 0) == Approx(4.0));
  CHECK(g(1, 1) == Approx(9.0));
  CHECK(g(2, 2) == Approx(16.0));
}

TEST_CASE("integer inverse of a unimodular rotation", "[math]") {
  Matrix3i rot; // 90 degrees about z
  rot << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  auto inv = math::integer_inverse(rot);
  REQUIRE(inv.has_value());
  CHECK((rot * *inv) == Matrix3i::Identity());
  // Non-unimodular -> no integer inverse.
  Matrix3i scale = Matrix3i::Identity() * 2;
  CHECK_FALSE(math::integer_inverse(scale).has_value());
}

TEST_CASE("similar matrix is b^-1 . a . b", "[math]") {
  Matrix3d a;
  a << 1, 2, 3, 0, 1, 4, 0, 0, 1;
  auto s = math::similar(a, Matrix3d::Identity(), 1e-10);
  REQUIRE(s.has_value());
  CHECK(s->isApprox(a, 1e-12));
}

TEST_CASE("frac displacement is the minimal image", "[math]") {
  Vector3d a(0.1, 0.9, 0.5);
  Vector3d b(0.2, 0.1, 0.5);
  Vector3d d = math::frac_displacement(a, b);
  CHECK(d[0] == Approx(0.1));
  CHECK(d[1] == Approx(0.2)); // 0.1 - 0.9 = -0.8 -> minimal image +0.2
  CHECK(d[2] == Approx(0.0));
}
