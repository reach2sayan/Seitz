// Oracle test for the layer-group dataset (spglib.c get_layer_dataset):
// the cppcrystal layer path must reproduce the reference spg_get_layer_dataset —
// the layer-group identity (negative hall number, layer-group number 1..80,
// symbol, point group), the symmetry operations, the per-atom Wyckoff /
// site-symmetry / equivalence data, and the standardized layer cell.
//
// The aperiodic axis (c in the standardized setting) is not periodic, so
// fractional comparisons of operation translations and standardized positions
// fold only the two periodic axes, never the aperiodic one.

#include "oracle.hpp"

#include <cppcrystal/dataset.hpp>

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

using cppcrystal::Cell;
using cppcrystal::Lattice;
using cppcrystal::Matrix3d;
using cppcrystal::Positions;
using cppcrystal::Result;
using cppcrystal::Types;
using cppcrystal::Vector3d;

Cell make_layer_cell(Matrix3d const &lattice,
                     std::vector<std::array<double, 3>> const &pos,
                     std::vector<int> const &types) {
  Positions p(static_cast<Eigen::Index>(pos.size()), 3);
  for (std::size_t i = 0; i < pos.size(); ++i) {
    p.row(static_cast<Eigen::Index>(i)) =
        Eigen::RowVector3d(pos[i][0], pos[i][1], pos[i][2]);
  }
  return Cell(Lattice{lattice}, p, types);
}

Matrix3d hexagonal_layer(double a, double c) {
  Matrix3d m;
  m.col(0) = Vector3d(a, 0, 0);
  m.col(1) = Vector3d(-a / 2, a * std::numbers::sqrt3 / 2, 0);
  m.col(2) = Vector3d(0, 0, c);
  return m;
}

Matrix3d orthorhombic_layer(double a, double b, double c) {
  Matrix3d m = Matrix3d::Zero();
  m(0, 0) = a;
  m(1, 1) = b;
  m(2, 2) = c;
  return m;
}

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
  Vector3d d = a - b;
  for (int j = 0; j < 3; ++j) {
    if (j != aperiodic_axis) {
      d[j] -= std::round(d[j]);
    }
  }
  return d.cwiseAbs().maxCoeff() <= symprec;
}

// Two standardized lattices describe the same cell up to a rigid rotation iff
// their metric tensors (L^T L) agree — spglib reports the orientation
// separately via std_rotation_matrix, and our layer idealization may orient the
// cell differently while reducing to the identical (a, b, c, angles).
bool same_metric(Matrix3d const &a, Matrix3d const &b) {
  Matrix3d const ga = a.transpose() * a;
  Matrix3d const gb = b.transpose() * b;
  return (ga - gb).cwiseAbs().maxCoeff() < 1e-3;
}

using cppcrystal::test::must;

void check(Cell const &cell, int aperiodic_axis, double symprec) {
  auto const got = must(cppcrystal::get_dataset(
      cell.with_periodicity(cppcrystal::aperiodic_along(aperiodic_axis)),
      {symprec}));
  auto const ref =
      cppcrystal::oracle::reference_layer_dataset(cell, aperiodic_axis, symprec);
  REQUIRE(ref.number != 0); // reference succeeded

  // Layer-group identity.
  CHECK(got.spacegroup_number == ref.number);
  CHECK(got.hall_number == ref.hall_number);
  CHECK(std::string(got.international_symbol) == ref.international);
  CHECK(std::string(got.hall_symbol) == ref.hall_symbol);
  CHECK(std::string(got.pointgroup_symbol) == ref.pointgroup);
  CHECK(got.aperiodic_axis.has_value());

  // Per input-atom Wyckoff / site-symmetry / equivalence data.
  REQUIRE(got.wyckoffs.size() == ref.wyckoffs.size());
  for (std::size_t i = 0; i < ref.wyckoffs.size(); ++i) {
    INFO("atom " << i);
    CHECK(got.wyckoffs[i] == ref.wyckoffs[i]);
    CHECK(got.site_symmetry_symbols[i] == ref.site_symmetry_symbols[i]);
    CHECK(got.equivalent_atoms[i] == ref.equivalent_atoms[i]);
  }

  // Symmetry operations, matched as a set (aperiodic-aware translation compare).
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
  CHECK(static_cast<int>(got.std_types.size()) == ref.n_std_atoms);
  CHECK(same_metric(got.std_lattice, ref.std_lattice));
  for (int r = 0; r < ref.n_std_atoms; ++r) {
    Vector3d const ref_pos = ref.std_positions.row(r).transpose();
    bool found = false;
    for (int g = 0; g < static_cast<int>(got.std_types.size()) && !found; ++g) {
      auto const ug = static_cast<std::size_t>(g);
      if (got.std_types[ug] != ref.std_types[static_cast<std::size_t>(r)]) {
        continue;
      }
      // Standardized positions are stored folded into [0, 1) on every axis; the
      // aperiodic axis is c = 2 in the standardized setting.
      found = overlap_layer(got.std_positions.row(g).transpose(), ref_pos,
                            /*aperiodic_axis=*/2, symprec);
    }
    INFO("standardized atom " << r << " not matched");
    CHECK(found);
  }
}

} // namespace

TEST_CASE("layer dataset: graphene p6/mmm", "[oracle][layer]") {
  Cell const cell = make_layer_cell(
      hexagonal_layer(2.46, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.0}, {2.0 / 3, 1.0 / 3, 0.0}}, {6, 6});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: graphene offset along c", "[oracle][layer]") {
  Cell const cell = make_layer_cell(
      hexagonal_layer(2.46, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.3}, {2.0 / 3, 1.0 / 3, 0.3}}, {6, 6});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: hexagonal BN p-6m2", "[oracle][layer]") {
  // Two distinct species on the honeycomb sites -> the mirror-only hexagonal
  // layer group.
  Cell const cell = make_layer_cell(
      hexagonal_layer(2.50, 15.0),
      {{1.0 / 3, 2.0 / 3, 0.0}, {2.0 / 3, 1.0 / 3, 0.0}}, {5, 7});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: square lattice p4/mmm", "[oracle][layer]") {
  Cell const cell =
      make_layer_cell(orthorhombic_layer(3.0, 3.0, 12.0), {{0.0, 0.0, 0.0}},
                      {1});
  check(cell, 2, 1e-5);
}

TEST_CASE("layer dataset: rectangular orthorhombic pmmm", "[oracle][layer]") {
  // The orthorhombic axis-choice loop is restricted to the two c-preserving
  // choices (abc, ba-c) with the layer axis-choice table.
  check(make_layer_cell(orthorhombic_layer(3.0, 5.0, 12.0), {{0.0, 0.0, 0.0}},
                        {1}),
        2, 1e-5);
}

TEST_CASE("layer dataset: centered rectangular (oc) lattice",
          "[oracle][layer]") {
  check(make_layer_cell(orthorhombic_layer(3.0, 5.0, 12.0),
                        {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}}, {1, 1}),
        2, 1e-5);
}

TEST_CASE("layer dataset: oblique monoclinic layers", "[oracle][layer]") {
  SECTION("single atom at the origin") {
    check(make_layer_cell(oblique_layer(3.0, 4.0, 105.0, 12.0),
                          {{0.0, 0.0, 0.0}}, {1}),
          2, 1e-5);
  }
  SECTION("single atom at a general position") {
    check(make_layer_cell(oblique_layer(3.0, 4.0, 105.0, 12.0),
                          {{0.2, 0.3, 0.0}}, {1}),
          2, 1e-5);
  }
  SECTION("two atoms related by a 2-fold") {
    check(make_layer_cell(oblique_layer(3.0, 4.0, 110.0, 12.0),
                          {{0.2, 0.1, 0.0}, {0.8, 0.9, 0.0}}, {1, 1}),
          2, 1e-5);
  }
}
