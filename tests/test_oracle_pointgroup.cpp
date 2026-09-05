#include "oracle.hpp"

#include "symmetry/pointgroup.hpp"
#include "symmetry/search.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace seitz;

namespace {
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

Cell hexagonal() {
  Matrix3d l;
  l.col(0) = Vector3d(3.0, 0.0, 0.0);
  l.col(1) = Vector3d(-1.5, 2.598076211, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, 5.0);
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Lattice{l}, pos, {0});
}
} // namespace

TEST_CASE("identify_point_group matches spg_get_pointgroup",
          "[oracle][pointgroup]") {
  for (Cell const &cell : {primitive_cubic(4.0), rutile(), hexagonal()}) {
    auto ops =
        symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}.operations();
    REQUIRE(ops);
    // Feed identical rotations to both implementations so axis selection and
    // de-duplication order match.
    auto const rotations = ops->rotations();
    auto ours = symmetry::identify_point_group<GroupFamily::space>(rotations);
    REQUIRE(ours);
    auto const ref = oracle::reference_pointgroup(rotations);

    INFO("ours = " << ours->pointgroup.symbol << " (" << ours->pointgroup.number
                   << "), ref = " << ref.symbol << " (" << ref.number << ")");
    CHECK(ours->pointgroup.number == ref.number);
    CHECK(std::string(ours->pointgroup.symbol) == ref.symbol);
    CHECK(ours->transformation == ref.transformation);
  }
}

TEST_CASE("crystal_class is aligned to the point-group numbering",
          "[pointgroup]") {
  // The documented invariant that consumers (e.g. switch on CrystalClass) rely
  // on: static_cast<CrystalClass>(number) round-trips through the table.
  for (int n = 1; n <= 32; ++n) {
    auto const pg = pointgroup_by_number(n);
    INFO("number = " << n << ", schoenflies = " << pg.schoenflies);
    CHECK(pg.crystal_class == static_cast<CrystalClass>(n));
    CHECK(pg.crystal_class != CrystalClass::none);
  }
  // Out-of-range numbers carry no class.
  CHECK(pointgroup_by_number(0).crystal_class == CrystalClass::none);
  CHECK(pointgroup_by_number(33).crystal_class == CrystalClass::none);
}
