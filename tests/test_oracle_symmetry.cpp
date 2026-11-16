#include "oracle.hpp"
#include <seitz/core/operation_set.hpp>

#include "symmetry/search.hpp"
#include "structures.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace seitz;

using seitz::test::primitive_cubic;
using seitz::test::rock_salt;
using seitz::test::rutile;
using seitz::test::same_operation_set;
using seitz::test::triclinic;

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
