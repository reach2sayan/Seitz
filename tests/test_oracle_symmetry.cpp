#include "oracle.hpp"
#include <cppcrystal/core/operation_set.hpp>

#include "symmetry/search.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;

namespace {
// Returns true iff the two operation collections are equal as sets (rotations
// exact, translations modulo the lattice within tolerance). spg_get_symmetry
// returns operations in a different order than we discover them.
bool same_operation_set(Operations const &a, Operations const &b, double tol) {
  if (a.size() != b.size())
    return false;
  for (auto const &oa : a) {
    bool matched = false;
    for (auto const &ob : b)
      if (same_operation(oa, ob, tol)) {
        matched = true;
        break;
      }
    if (!matched)
      return false;
  }
  return true;
}

Cell primitive_cubic(double a) {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Lattice{Matrix3d::Identity() * a}, pos, {0});
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
  return Cell(Lattice{l}, pos, {0, 0, 1, 1, 1, 1});
}

Cell rock_salt(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(8, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  pos.row(2) << 0.5, 0.0, 0.5;
  pos.row(3) << 0.0, 0.5, 0.5;
  pos.row(4) << 0.5, 0.5, 0.5;
  pos.row(5) << 0.0, 0.0, 0.5;
  pos.row(6) << 0.0, 0.5, 0.0;
  pos.row(7) << 0.5, 0.0, 0.0;
  return Cell(Lattice{lattice}, pos, {0, 0, 0, 0, 1, 1, 1, 1});
}

// Low-symmetry P1-ish cell with two unlike atoms in general positions.
Cell triclinic() {
  Matrix3d l;
  l.col(0) = Vector3d(4.0, 0.0, 0.0);
  l.col(1) = Vector3d(0.7, 4.3, 0.0);
  l.col(2) = Vector3d(0.5, 0.9, 5.1);
  Positions pos(2, 3);
  pos.row(0) << 0.1, 0.2, 0.3;
  pos.row(1) << 0.6, 0.55, 0.27;
  return Cell(Lattice{l}, pos, {0, 1});
}
} // namespace

TEST_CASE("find_symmetry matches spg_get_symmetry as a set",
          "[oracle][symmetry]") {
  for (Cell const &cell :
       {primitive_cubic(4.0), rutile(), rock_salt(5.6), triclinic()}) {
    auto ours =
        symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}.operations();
    REQUIRE(ours);
    auto const ref = oracle::reference_symmetry(cell, 1e-5);
    INFO("ours = " << ours->size() << ", reference = " << ref.size());
    CHECK(same_operation_set(*ours, ref, 1e-4));
  }
}
