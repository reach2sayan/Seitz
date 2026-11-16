// Oracle test for the symmetry -> space-group pipeline (the foundation the
// magnetic determination builds on): OperationSet::spacegroup
// (= prm_get_primitive_symmetry + spa_search_spacegroup_with_symmetry) must
// reproduce spg_get_spacegroup_type_from_symmetry. The same operation set is
// fed to both sides, so this tests the pipeline, not the symmetry search.

#include "oracle.hpp"
#include "structures.hpp"
#include <seitz/core/operation_set.hpp>

#include "spacegroup/spacegroup.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using seitz::test::make_cell;
using seitz::test::cubic;
using seitz::test::hexagonal;
using seitz::test::tetragonal;

using seitz::Cell;
using seitz::Lattice;
using seitz::Matrix3d;
using seitz::Operations;
using seitz::Positions;
using seitz::oracle::reference_symmetry;

// Compare port OperationSet::spacegroup<conventional> against the reference
// spg_get_spacegroup_type_from_symmetry, both fed the same operation set.
void check_conventional(Cell const &cell, double symprec) {
  auto const ops = reference_symmetry(cell, symprec);
  REQUIRE(!ops.empty());

  auto const got = ops.spacegroup(cell.lattice(), {symprec});
  REQUIRE(got);

  seitz::oracle::COperations c_ops(ops);
  double lat[3][3];
  seitz::oracle::to_c_lattice(lat, cell.lattice().matrix());
  SpglibSpacegroupType const ref = spg_get_spacegroup_type_from_symmetry(
      c_ops.rotations(), c_ops.translations(), c_ops.size(), lat, symprec);

  REQUIRE(ref.number != 0);
  REQUIRE(got->type().number == ref.number);
  REQUIRE(got->hall.index() == ref.hall_number);
  REQUIRE(got->type().international_short ==
          std::string(ref.international_short));
}

} // namespace

TEST_CASE("OperationSet::spacegroup matches reference (conventional)",
          "[oracle][spacegroup][from_symmetry]") {
  double const s = 1e-5;
  SECTION("Pm-3m (primitive cubic)") {
    check_conventional(make_cell(cubic(4.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("Im-3m (body-centred cubic)") {
    check_conventional(
        make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0}), s);
  }
  SECTION("Fm-3m (face-centred cubic)") {
    check_conventional(
        make_cell(cubic(4.0),
                  {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}},
                  {0, 0, 0, 0}),
        s);
  }
  SECTION("P6/mmm (hexagonal)") {
    check_conventional(make_cell(hexagonal(3.0, 5.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("P4/mmm (tetragonal)") {
    check_conventional(make_cell(tetragonal(3.0, 5.0), {{0, 0, 0}}, {0}), s);
  }
  SECTION("rutile-like (tetragonal, 2 species)") {
    check_conventional(make_cell(tetragonal(4.6, 3.0),
                                 {{0, 0, 0},
                                  {0.5, 0.5, 0.5},
                                  {0.3, 0.3, 0.0},
                                  {0.7, 0.7, 0.0},
                                  {0.8, 0.2, 0.5},
                                  {0.2, 0.8, 0.5}},
                                 {0, 0, 1, 1, 1, 1}),
                       s);
  }
}

// Exercises the LatticeSetting::primitive branch (identity primitive lattice),
// mirroring spg_get_hall_number_from_symmetry.
TEST_CASE("OperationSet::spacegroup matches reference (primitive)",
          "[oracle][spacegroup][from_symmetry]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0});
  auto const ops = reference_symmetry(cell, s);
  REQUIRE(!ops.empty());

  seitz::oracle::COperations c_ops(ops);
  int const ref_hall = spg_get_hall_number_from_symmetry(
      c_ops.rotations(), c_ops.translations(), c_ops.size(), s);

  auto const got = ops.spacegroup<seitz::LatticeSetting::primitive>(
      *Lattice::from_basis(Matrix3d::Identity()), {s});

  REQUIRE((ref_hall == 0) == (!got.has_value()));
  if (got) {
    REQUIRE(got->hall.index() == ref_hall);
  }
}
