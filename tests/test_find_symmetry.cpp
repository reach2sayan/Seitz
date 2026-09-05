#include "symmetry/search.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace seitz;

namespace {
Cell primitive_cubic(double a) {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Lattice{Matrix3d::Identity() * a}, pos, {0});
}

Cell rock_salt(double a) {
  // Conventional NaCl cell: 4 Na + 4 Cl (Fm-3m). Non-primitive -> many ops.
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
} // namespace

TEST_CASE("primitive cubic has the full 48-operation point group",
          "[symmetry]") {
  Cell const cell = primitive_cubic(4.0);
  auto ops =
      symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}.operations();
  REQUIRE(ops);
  CHECK(ops->size() == 48);
  // All rotations unimodular; identity present with zero translation.
  bool has_identity = false;
  for (auto const &op : *ops) {
    CHECK(std::abs(op.rotation.determinant()) == 1);
    if (op.is_identity_rotation() &&
        op.translation.cwiseAbs().maxCoeff() < 1e-9)
      has_identity = true;
  }
  CHECK(has_identity);
}

TEST_CASE("conventional rock-salt cell has 192 operations", "[symmetry]") {
  // Fm-3m: 48 point ops x 4 centering translations.
  Cell const cell = rock_salt(5.6);
  auto ops =
      symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}.operations();
  REQUIRE(ops);
  CHECK(ops->size() == 192);
}

TEST_CASE("lattice point group of a tetragonal cell has 16 operations",
          "[symmetry]") {
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.0;
  l(1, 1) = 4.0;
  l(2, 2) = 6.0;
  Cell const cell(Lattice{l}, Positions(0, 3), {});
  auto ps = symmetry::SymmetrySearch<GroupFamily::space>{
      cell,
      {1e-5}}.lattice_symmetry();
  REQUIRE(ps);
  CHECK(ps->size() == 16); // 4/mmm
}
