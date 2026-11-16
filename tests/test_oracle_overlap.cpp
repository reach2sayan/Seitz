#include "oracle.hpp"

#include "core/overlap.hpp"
#include "structures.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace seitz;

using seitz::test::bcc;
using seitz::test::rutile;

TEST_CASE("every reference symmetry op overlaps under our checker",
          "[oracle][overlap]") {
  for (Cell const &cell : {bcc(3.1), rutile()}) {
    OverlapChecker checker(cell, 1e-5);
    auto const ops = oracle::reference_symmetry(cell, 1e-5);
    REQUIRE(!ops.empty());
    for (auto const &op : ops) {
      INFO("rotation det = " << op.rotation.determinant());
      CHECK(checker.check_total_overlap(op.translation, op.rotation));
    }
  }
}

TEST_CASE("a corrupted reference op is rejected", "[oracle][overlap]") {
  Cell const cell = rutile();
  OverlapChecker checker(cell, 1e-5);
  auto ops = oracle::reference_symmetry(cell, 1e-5);
  REQUIRE(ops.size() > 1);
  // Shift a real operation's translation off by a clearly non-symmetric amount.
  SymmetryOperation bad = ops[0];
  bad.translation += Vector3d(0.13, 0.0, 0.0);
  CHECK_FALSE(checker.check_total_overlap(bad.translation, bad.rotation));
}
