// Oracle test for the irreducible / stabilized reciprocal mesh (kpoint.c):
// ir_reciprocal_mesh, stabilized_reciprocal_mesh, and grid_points_by_rotations
// must reproduce spg_get_dense_ir_reciprocal_mesh /
// spg_get_dense_stabilized_reciprocal_mesh / spg_get_dense_grid_points_by_
// rotations.

#include "oracle.hpp"
#include "structures.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/kpoint/mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace {

using seitz::test::make_cell;
using seitz::test::cubic;
using seitz::test::hexagonal;

using seitz::Cell;
using seitz::Lattice;
using seitz::Matrix3d;
using seitz::Matrix3i;
using seitz::Positions;
using seitz::Vector3d;
using seitz::Vector3i;

seitz::kpoint::Mesh sampling(Vector3i const &divisions,
                                  Vector3i const &is_shift) {
  auto const mesh = seitz::kpoint::Mesh::of(
      {divisions[0], divisions[1], divisions[2]},
      {is_shift[0] != 0, is_shift[1] != 0, is_shift[2] != 0});
  REQUIRE(mesh);
  return *mesh;
}

void check_ir(Cell const &cell, Vector3i const &divisions,
              Vector3i const &is_shift, bool time_reversal, double symprec) {
  auto const analyzer =
      seitz::analysis::SymmetryAnalyzer::from_cell(cell, {symprec});
  auto const got =
      analyzer.reciprocal_mesh(sampling(divisions, is_shift),
                               time_reversal ? seitz::TimeReversal::on
                                             : seitz::TimeReversal::off);
  REQUIRE(got);
  auto const ref = seitz::oracle::reference_ir_reciprocal_mesh(
      cell, divisions, is_shift, time_reversal, symprec);

  CHECK(got->num_irreducible() == ref.num_ir);
  REQUIRE(got->mesh().size() == ref.grid_address.size());
  REQUIRE(got->mapping().size() == ref.ir_mapping_table.size());
  for (auto const [i, point] :
       std::views::zip(got->mesh().addresses(), got->mapping(),
                       ref.grid_address, ref.ir_mapping_table) |
           std::views::enumerate) {
    auto const &[address, mapped, ref_address, ref_mapped] = point;
    INFO("grid point " << static_cast<std::size_t>(i));
    CHECK(Vector3i(address[0], address[1], address[2]) == ref_address);
    CHECK(mapped == ref_mapped);
  }
}

} // namespace

TEST_CASE("ir reciprocal mesh matches reference (normal path)",
          "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const pm3m = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  Cell const im3m = make_cell(cubic(4.0), {{0, 0, 0}, {0.5, 0.5, 0.5}}, {0, 0});

  SECTION("gamma-centred") {
    check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), true, s);
  }
  SECTION("shifted") {
    check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(1, 1, 1), true, s);
  }
  SECTION("odd mesh") {
    check_ir(pm3m, Vector3i(3, 3, 3), Vector3i(0, 0, 0), true, s);
  }
  SECTION("no time reversal") {
    check_ir(pm3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), false, s);
  }
  SECTION("bcc") {
    check_ir(im3m, Vector3i(4, 4, 4), Vector3i(0, 0, 0), true, s);
  }
}

TEST_CASE("ir reciprocal mesh matches reference (distortion path)",
          "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const p6mmm = make_cell(hexagonal(3.0, 5.0), {{0, 0, 0}}, {0});
  SECTION("gamma-centred") {
    check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 0), true, s);
  }
  SECTION("shifted c") {
    check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 1), true, s);
  }
  SECTION("no time reversal") {
    check_ir(p6mmm, Vector3i(6, 6, 4), Vector3i(0, 0, 0), false, s);
  }
}

TEST_CASE("stabilized reciprocal mesh matches reference", "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  auto const rotations =
      seitz::oracle::reference_symmetry(cell, s).rotations();
  Vector3i const mesh(4, 4, 4);
  Vector3i const shift(0, 0, 0);

  auto compare = [&](std::vector<Vector3d> const &qpoints, bool tr) {
    auto const got =
        seitz::kpoint::ReciprocalMesh::from_rotations(
            sampling(mesh, shift), rotations,
            tr ? seitz::TimeReversal::on : seitz::TimeReversal::off)
            .stabilized(qpoints);
    auto const ref = seitz::oracle::reference_stabilized_reciprocal_mesh(
        mesh, shift, tr, rotations, qpoints);
    CHECK(got.num_irreducible() == ref.num_ir);
    REQUIRE(got.mapping().size() == ref.ir_mapping_table.size());
    CHECK(std::ranges::equal(got.mapping(), ref.ir_mapping_table));
  };

  SECTION("gamma only (full symmetry)") {
    compare({Vector3d(0.0, 0.0, 0.0)}, true);
  }
  SECTION("q breaking symmetry") { compare({Vector3d(0.5, 0.0, 0.0)}, true); }
  SECTION("q set, no time reversal") {
    compare({Vector3d(0.5, 0.0, 0.0), Vector3d(0.0, 0.5, 0.0)}, false);
  }
}

TEST_CASE("grid points by rotations match reference", "[oracle][kpoint]") {
  double const s = 1e-5;
  Cell const cell = make_cell(cubic(4.0), {{0, 0, 0}}, {0});
  auto const rotations =
      seitz::oracle::reference_symmetry(cell, s).rotations();
  Vector3i const mesh(4, 4, 4);

  for (Vector3i const &shift : {Vector3i(0, 0, 0), Vector3i(1, 1, 1)}) {
    auto const reciprocal = seitz::kpoint::ReciprocalMesh::from_rotations(
        sampling(mesh, shift), rotations, seitz::TimeReversal::off);
    std::vector<Matrix3i> const rot_reciprocal(reciprocal.rotations().begin(),
                                               reciprocal.rotations().end());
    for (Vector3i const &address : {Vector3i(1, 2, 0), Vector3i(2, 0, 1)}) {
      INFO("shift " << shift.transpose() << " address " << address.transpose());
      auto const got =
          reciprocal.images_of({address[0], address[1], address[2]});
      auto const ref = seitz::oracle::reference_grid_points_by_rotations(
          address, rot_reciprocal, mesh, shift);
      CHECK(got == ref);
    }
  }
}
