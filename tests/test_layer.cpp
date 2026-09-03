// Layer-group tests (no reference-spglib oracle needed): known 2D-periodic
// structures with hand-checked layer-group numbers, symbols and Wyckoff data.
// The oracle comparison against spg_get_layer_dataset lives in the gated
// oracle suite; these cases pin the behavior in the default build.
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/data/spg_database.hpp>

#include "helpers.hpp"

#include <boost/leaf.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>
#include <utility>

using namespace cppcrystal;

namespace {
using cppcrystal::test::layer_hall;
using cppcrystal::test::space_hall;
using cppcrystal::test::must;

// A hexagonal in-plane lattice (a=b, gamma=120) with a large vacuum gap along c,
// the aperiodic axis. Columns are basis vectors.
Matrix3d hexagonal_layer(double a, double c) {
  Matrix3d m;
  m.col(0) = Vector3d(a, 0, 0);
  m.col(1) = Vector3d(-a / 2, a * std::numbers::sqrt3 / 2, 0);
  m.col(2) = Vector3d(0, 0, c);
  return m;
}
} // namespace

TEST_CASE("layer database: negative-hall metadata + ops", "[layer][database]") {
  // constexpr metadata works at compile time for layer settings.
  STATIC_REQUIRE(data::spacegroup_type(layer_hall(1)).number == 1);
  STATIC_REQUIRE(data::spacegroup_type(layer_hall(116)).number == 80);
  STATIC_REQUIRE(data::spacegroup_type(layer_hall(116)).pointgroup_number == 27);
  STATIC_REQUIRE(data::spacegroup_type(space_hall(1)).number == 1); // 3D path intact

  // Layer operations are capped at 24 (vs 48 for 3D point groups).
  REQUIRE(data::operations_from_database(layer_hall(116)).size() == 24); // p6/mmm layer
  REQUIRE(data::operations_from_database(layer_hall(1)).size() == 1);     // p1 layer
}

TEST_CASE("graphene is layer group p6/mmm (LG 80)", "[layer]") {
  Matrix3d const lat = hexagonal_layer(2.46, 15.0);
  Positions pos(2, 3);
  pos.row(0) << 1.0 / 3, 2.0 / 3, 0.0; // honeycomb sites in the z=0 plane
  pos.row(1) << 2.0 / 3, 1.0 / 3, 0.0;
  Cell const cell{Lattice{lat}, pos, Types{6, 6}};

  auto ds = must(test::dataset_of(cell.with_periodicity(aperiodic_along(2)), {1e-4}));
  REQUIRE(data::spacegroup_type(ds.hall).number == 80);  // layer-group number
  REQUIRE(ds.hall == layer_hall(116));      // negative-hall convention
  REQUIRE(data::spacegroup_type(ds.hall).international_short == "p6/mmm");
  REQUIRE(aperiodic_axis(ds.standardized.periodicity()) == 2);
  REQUIRE(ds.operations.size() == 24);
  // Both carbons sit on the same Wyckoff orbit with -6m2 site symmetry.
  REQUIRE(ds.sites.size() == 2);
  REQUIRE(ds.sites[0].wyckoff == ds.sites[1].wyckoff);
  REQUIRE(ds.sites[0].site_symmetry == "-6m2");
}

TEST_CASE("graphene offset along the aperiodic axis still resolves", "[layer]") {
  // The layer sits at an arbitrary z (not the database's z=0 plane). The
  // origin-shift aperiodic (c) component must re-center it; the result is
  // independent of the offset.
  Matrix3d const lat = hexagonal_layer(2.46, 15.0);
  for (double z : {0.5, 0.27}) {
    Positions pos(2, 3);
    pos.row(0) << 1.0 / 3, 2.0 / 3, z;
    pos.row(1) << 2.0 / 3, 1.0 / 3, z;
    Cell const cell{Lattice{lat}, pos, Types{6, 6}};
    auto ds = must(test::dataset_of(cell.with_periodicity(aperiodic_along(2)), {1e-4}));
    INFO("z offset = " << z);
    REQUIRE(data::spacegroup_type(ds.hall).number == 80);
    REQUIRE(data::spacegroup_type(ds.hall).international_short == "p6/mmm");
    REQUIRE(ds.sites[0].site_symmetry == "-6m2");
  }
}

TEST_CASE("square lattice single atom is p4/mmm (LG 61)", "[layer]") {
  Matrix3d lat;
  lat.col(0) = Vector3d(1, 0, 0);
  lat.col(1) = Vector3d(0, 1, 0);
  lat.col(2) = Vector3d(0, 0, 8);
  Positions pos(1, 3);
  pos.row(0) << 0, 0, 0;
  Cell const cell{Lattice{lat}, pos, Types{1}};

  auto ds = must(test::dataset_of(cell.with_periodicity(aperiodic_along(2)), {1e-4}));
  REQUIRE(data::spacegroup_type(ds.hall).number == 61);
  REQUIRE(data::spacegroup_type(ds.hall).international_short == "p4/mmm");
  REQUIRE(ds.sites[0].site_symmetry == "4/mmm");
}

TEST_CASE("SymmetryAnalyzer auto-routes layer cells", "[layer][analysis]") {
  Matrix3d const lat = hexagonal_layer(2.46, 15.0);
  Positions pos(2, 3);
  pos.row(0) << 1.0 / 3, 2.0 / 3, 0.0;
  pos.row(1) << 2.0 / 3, 1.0 / 3, 0.0;
  Cell const cell{Lattice{lat}, pos, Types{6, 6}};

  auto analyzer = analysis::SymmetryAnalyzer::from_cell(
      cell.with_periodicity(aperiodic_along(2)), {1e-4});
  REQUIRE(data::spacegroup_type(must(analyzer.hall())).number == 80);
  REQUIRE(must(analyzer.hall()) == layer_hall(116));
}

TEST_CASE("a 3D cell is unaffected by the layer code paths", "[layer]") {
  // Simple cubic single atom -> Pm-3m (221); the aperiodic branches must not
  // perturb the 3D result.
  Positions pos(1, 3);
  pos.row(0) << 0, 0, 0;
  Cell const cell{Lattice{Matrix3d::Identity()}, pos, Types{1}};
  auto ds = must(test::dataset_of(cell));
  REQUIRE(data::spacegroup_type(ds.hall).number == 221);
  REQUIRE_FALSE(aperiodic_axis(ds.standardized.periodicity()).has_value());
}
