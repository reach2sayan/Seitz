#include "oracle.hpp"

#include <cppcrystal/reduce/delaunay.hpp>
#include <cppcrystal/reduce/niggli.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>

using namespace cppcrystal;

namespace {
// A spread of lattices: cubic, a primitive-rhombohedral-ish cell, a sheared
// monoclinic, and a strongly oblique triclinic. Columns = basis vectors.
std::array<Matrix3d, 4> test_lattices() {
  std::array<Matrix3d, 4> v;
  v[0] = Matrix3d::Identity() * 3.7;
  v[1].col(0) = Vector3d(2.0, 0.0, 0.0);
  v[1].col(1) = Vector3d(1.0, 2.0, 0.0);
  v[1].col(2) = Vector3d(1.0, 1.0, 2.0);
  v[2].col(0) = Vector3d(4.0, 0.0, 0.0);
  v[2].col(1) = Vector3d(0.0, 5.0, 0.0);
  v[2].col(2) = Vector3d(1.3, 0.0, 6.0); // monoclinic beta
  v[3].col(0) = Vector3d(3.1, 0.2, 0.1);
  v[3].col(1) = Vector3d(1.7, 3.4, 0.3);
  v[3].col(2) = Vector3d(0.9, 1.1, 3.9);
  return v;
}

// The reduced lattice is unique up to nothing here (both algorithms apply
// integer column ops to the SAME input), so we expect the Cartesian vectors to
// match.
void check_close(Matrix3d const &a, Matrix3d const &b) {
  CHECK((a - b).cwiseAbs().maxCoeff() < 1e-8);
}
} // namespace

TEST_CASE("delaunay reduction matches the reference", "[oracle][reduce]") {
  for (Matrix3d const &lat : test_lattices()) {
    auto ours = reduce::delaunay_reduce(lat, 1e-5);
    REQUIRE(ours);
    check_close(*ours, oracle::reference_delaunay(lat, 1e-5));
  }
}

TEST_CASE("niggli reduction matches the reference", "[oracle][reduce]") {
  for (Matrix3d const &lat : test_lattices()) {
    auto ours = reduce::niggli_reduce(lat, 1e-5);
    REQUIRE(ours);
    check_close(*ours, oracle::reference_niggli(lat, 1e-5));
  }
}
