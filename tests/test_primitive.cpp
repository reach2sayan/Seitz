#include "symmetry/primitive.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace seitz;
using Catch::Approx;

namespace {
Cell bcc_conventional(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0});
}

Cell fcc_conventional(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(4, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  pos.row(2) << 0.5, 0.0, 0.5;
  pos.row(3) << 0.0, 0.5, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0, 0, 0});
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
} // namespace

TEST_CASE("bcc conventional cell reduces to a 1-atom primitive cell",
          "[primitive]") {
  Cell const cell = bcc_conventional(3.0);
  auto p = symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
  REQUIRE(p);
  CHECK(p->cell.size() == 1);
  CHECK(p->cell.lattice().volume() ==
        Approx(27.0 / 2.0)); // half the conventional volume
}

TEST_CASE("fcc conventional cell reduces to a 1-atom primitive cell",
          "[primitive]") {
  Cell const cell = fcc_conventional(4.0);
  auto p = symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
  REQUIRE(p);
  CHECK(p->cell.size() == 1);
  CHECK(p->cell.lattice().volume() == Approx(64.0 / 4.0));
}

TEST_CASE("rock-salt conventional cell reduces to a 2-atom primitive cell",
          "[primitive]") {
  Cell const cell = rock_salt(5.6);
  auto p = symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
  REQUIRE(p);
  CHECK(p->cell.size() == 2); // one of each species
  // mapping_table covers all 8 input atoms, into 2 primitive atoms.
  CHECK(p->mapping_table.size() == 8);
}

TEST_CASE("an already-primitive cell is returned with one atom",
          "[primitive]") {
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  Cell const cell(Lattice{Matrix3d::Identity() * 4.0}, pos, {0});
  auto p = symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
  REQUIRE(p);
  CHECK(p->cell.size() == 1);
}
