#include <spglib/symmetry/find_symmetry.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace spglib;

namespace {
Cell primitive_cubic(double a) {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Matrix3d::Identity() * a, pos, {0});
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
  return Cell(lattice, pos, {0, 0, 0, 0, 1, 1, 1, 1});
}
} // namespace

TEST_CASE("primitive cubic has the full 48-operation point group",
          "[symmetry]") {
  auto ops = symmetry::find_symmetry(primitive_cubic(4.0), 1e-5);
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
  auto ops = symmetry::find_symmetry(rock_salt(5.6), 1e-5);
  REQUIRE(ops);
  CHECK(ops->size() == 192);
}

TEST_CASE("lattice point group of a tetragonal cell has 16 operations",
          "[symmetry]") {
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.0;
  l(1, 1) = 4.0;
  l(2, 2) = 6.0;
  auto ps = symmetry::lattice_symmetry(Cell(l, Positions(0, 3), {}), 1e-5);
  REQUIRE(ps);
  CHECK(ps->size() == 16); // 4/mmm
}
