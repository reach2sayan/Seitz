#include "oracle.hpp"

#include "symmetry/primitive.hpp"
#include <cppcrystal/core/lattice.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using Catch::Approx;

namespace {
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

// Two formula units of a tetragonal body-centered structure (I-centered).
Cell body_centered_tetragonal() {
  Matrix3d lattice = Matrix3d::Zero();
  lattice(0, 0) = 4.0;
  lattice(1, 1) = 4.0;
  lattice(2, 2) = 6.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0});
}
} // namespace

TEST_CASE("find_primitive matches spg_find_primitive (size, lattice, volume)",
          "[oracle][primitive]") {
  for (Cell const &cell :
       {fcc_conventional(4.0), rock_salt(5.6), body_centered_tetragonal()}) {
    auto ours =
        symmetry::PrimitiveFinder<GroupFamily::space>{cell, {1e-5}}.find();
    REQUIRE(ours);
    Cell const ref = oracle::reference_find_primitive(cell, 1e-5);

    INFO("ours size = " << ours->cell.size() << ", ref size = " << ref.size());
    REQUIRE(ours->cell.size() == ref.size());
    CHECK(ours->cell.lattice().volume() == Approx(ref.lattice().volume()));
    // Both describe the same lattice; compare the Niggli-reduced metric tensor,
    // which is invariant to basis choice and rigid rotation.
    auto const gours = ours->cell.lattice().niggli(1e-5);
    auto const gref = ref.lattice().niggli(1e-5);
    REQUIRE(gours);
    REQUIRE(gref);
    Matrix3d const mours = gours->metric();
    Matrix3d const mref = gref->metric();
    CHECK((mours - mref).cwiseAbs().maxCoeff() < 1e-6);
  }
}
