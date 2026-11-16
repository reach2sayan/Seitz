#include "oracle.hpp"

#include "symmetry/primitive.hpp"
#include "structures.hpp"
#include <seitz/core/lattice.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace seitz;
using Catch::Approx;

using seitz::test::bcc_tetragonal;
using seitz::test::fcc;
using seitz::test::rock_salt;

TEST_CASE("find_primitive matches spg_find_primitive (size, lattice, volume)",
          "[oracle][primitive]") {
  for (Cell const &cell :
       {fcc(4.0), rock_salt(5.6), bcc_tetragonal()}) {
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
