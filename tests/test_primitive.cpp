#include "symmetry/primitive.hpp"
#include "structures.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace seitz;
using Catch::Approx;

using seitz::test::bcc;
using seitz::test::fcc;
using seitz::test::rock_salt;

TEST_CASE("bcc conventional cell reduces to a 1-atom primitive cell",
          "[primitive]") {
  Cell const cell = bcc(3.0);
  auto p = symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
  REQUIRE(p);
  CHECK(p->cell.size() == 1);
  CHECK(p->cell.lattice().volume() ==
        Approx(27.0 / 2.0)); // half the conventional volume
}

TEST_CASE("fcc conventional cell reduces to a 1-atom primitive cell",
          "[primitive]") {
  Cell const cell = fcc(4.0);
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
