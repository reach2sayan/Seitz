// Oracle test for the layer-group dataset (spglib.c get_layer_dataset):
// the seitz layer path must reproduce the reference spg_get_layer_dataset
// — the layer-group identity (negative hall number, layer-group number 1..80,
// symbol, point group), the symmetry operations, the per-atom Wyckoff /
// site-symmetry / equivalence data, and the standardized layer cell.
//
// The aperiodic axis (c in the standardized setting) is not periodic, so
// fractional comparisons of operation translations and standardized positions
// fold only the two periodic axes, never the aperiodic one.

#include "oracle.hpp"
#include "structures.hpp"

#include <seitz/analysis/symmetry_analyzer.hpp>

#include <seitz/data/spg_database.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using seitz::test::make_cell;
using seitz::test::hexagonal;
using seitz::test::orthorhombic;

using seitz::Cell;
using seitz::Lattice;
using seitz::Matrix3d;
using seitz::Positions;
using seitz::Result;
using seitz::Types;
using seitz::Vector3d;
using seitz::test::layer_hall;

// Oblique 2D lattice (a != b, gamma != 90) with a vacuum gap along the
// perpendicular aperiodic axis c — the monoclinic/triclinic layer setting.
Matrix3d oblique_layer(double a, double b, double gamma_deg, double c) {
  double const g = gamma_deg * std::numbers::pi / 180.0;
  Matrix3d m;
  m.col(0) = Vector3d(a, 0, 0);
  m.col(1) = Vector3d(b * std::cos(g), b * std::sin(g), 0);
  m.col(2) = Vector3d(0, 0, c);
  return m;
}

// Fold a fractional displacement on the two periodic axes only (the aperiodic
// axis is not periodic).
bool overlap_layer(Vector3d const &a, Vector3d const &b, int aperiodic_axis,
                   double symprec) {
  return seitz::minimal_image(Vector3d(a - b),
                              seitz::aperiodic_along(aperiodic_axis))
             .cwiseAbs()
             .maxCoeff() <= symprec;
}

using seitz::test::must;

void check(Cell const &cell, int aperiodic_axis, double symprec) {
  auto const got = must(seitz::test::dataset_of(
      cell.with_periodicity(seitz::aperiodic_along(aperiodic_axis)),
      {symprec}));
  auto const ref = seitz::oracle::reference_layer_dataset(
      cell, aperiodic_axis, symprec);
  REQUIRE(ref.number != 0); // reference succeeded

  // Layer-group identity.
  CHECK(seitz::data::spacegroup_type(got.hall).number == ref.number);
  CHECK(got.hall == layer_hall(-ref.hall_number));
  CHECK(std::string(
            seitz::data::spacegroup_type(got.hall).international_short) ==
        ref.international);
  CHECK(std::string(seitz::data::spacegroup_type(got.hall).hall_symbol) ==
        ref.hall_symbol);
  CHECK(std::string(
            seitz::pointgroup_by_number(
                seitz::data::spacegroup_type(got.hall).pointgroup_number)
                .symbol) == ref.pointgroup);
  CHECK(seitz::aperiodic_axis(got.standardized.periodicity()).has_value());

  // Per input-atom Wyckoff / site-symmetry / equivalence data.
  seitz::oracle::check_sites(got.sites, ref);

  // Symmetry operations, matched as a set (aperiodic-aware translation
  // compare).
  REQUIRE(static_cast<int>(got.operations.size()) == ref.n_operations);
  for (auto const &rop : ref.operations) {
    bool found = false;
    for (auto const &gop : got.operations) {
      if (gop.rotation == rop.rotation &&
          overlap_layer(gop.translation, rop.translation, aperiodic_axis,
                        symprec)) {
        found = true;
        break;
      }
    }
    INFO("reference operation not found in our result");
    CHECK(found);
  }

  // Standardized layer cell: lattice + atoms (matched as a set), aperiodic c.
  CHECK(static_cast<int>(got.standardized.types().size()) == ref.n_std_atoms);
  // spglib reports the orientation separately via std_rotation_matrix, and our
  // layer idealization may orient the cell differently while reducing to the
  // identical (a, b, c, angles): compare metrics.
  CHECK(seitz::test::same_metric(got.standardized.lattice().matrix(),
                                 ref.std_lattice, 1e-3));
  for (int r = 0; r < ref.n_std_atoms; ++r) {
    Vector3d const ref_pos = ref.std_positions.row(r).transpose();
    bool found = false;
    for (int g = 0;
         g < static_cast<int>(got.standardized.types().size()) && !found; ++g) {
      auto const ug = static_cast<std::size_t>(g);
      if (got.standardized.types()[ug] !=
          ref.std_types[static_cast<std::size_t>(r)]) {
        continue;
      }
      // Standardized positions are stored folded into [0, 1) on every axis; the
      // aperiodic axis is c = 2 in the standardized setting.
      found = overlap_layer(got.standardized.positions().row(g).transpose(),
                            ref_pos,
                            /*aperiodic_axis=*/2, symprec);
    }
    INFO("standardized atom " << r << " not matched");
    CHECK(found);
  }
}

} // namespace

TEST_CASE("layer dataset: graphene p6/mmm", "[oracle][layer]") {
  Cell const cell = make_cell(
      hexagonal(2.46, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.0}, {2.0 / 3, 1.0 / 3, 0.0}}, {6, 6});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: graphene offset along c", "[oracle][layer]") {
  Cell const cell = make_cell(
      hexagonal(2.46, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.3}, {2.0 / 3, 1.0 / 3, 0.3}}, {6, 6});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: hexagonal BN p-6m2", "[oracle][layer]") {
  // Two distinct species on the honeycomb sites -> the mirror-only hexagonal
  // layer group.
  Cell const cell = make_cell(
      hexagonal(2.50, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.0}, {2.0 / 3, 1.0 / 3, 0.0}}, {5, 7});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: square lattice p4/mmm", "[oracle][layer]") {
  Cell const cell = make_cell(orthorhombic(3.0, 3.0, 12.0),
                                    {{0.0, 0.0, 0.0}}, {1});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: rectangular orthorhombic pmmm", "[oracle][layer]") {
  // The orthorhombic axis-choice loop is restricted to the two c-preserving
  // choices (abc, ba-c) with the layer axis-choice table.
  check(make_cell(orthorhombic(3.0, 5.0, 12.0), {{0.0, 0.0, 0.0}},
                        {1}),
        2, 1e-5);
}

TEST_CASE("layer dataset: centered rectangular (oc) lattice",
          "[oracle][layer]") {
  check(make_cell(orthorhombic(3.0, 5.0, 12.0),
                        {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}}, {1, 1}),
        2, 1e-5);
}

TEST_CASE("layer dataset: oblique monoclinic layers", "[oracle][layer]") {
  SECTION("single atom at the origin") {
    check(make_cell(oblique_layer(3.0, 4.0, 105.0, 12.0),
                          {{0.0, 0.0, 0.0}}, {1}),
          2, 1e-5);
  }
  SECTION("single atom at a general position") {
    check(make_cell(oblique_layer(3.0, 4.0, 105.0, 12.0),
                          {{0.2, 0.3, 0.0}}, {1}),
          2, 1e-5);
  }
  SECTION("two atoms related by a 2-fold") {
    check(make_cell(oblique_layer(3.0, 4.0, 110.0, 12.0),
                          {{0.2, 0.1, 0.0}, {0.8, 0.9, 0.0}}, {1, 1}),
          2, 1e-5);
  }
}
