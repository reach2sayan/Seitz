#include "oracle.hpp"

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/dataset.hpp>
#include <cppcrystal/refine/refinement.hpp>
#include <cppcrystal/refine/standardize.hpp>
#include <cppcrystal/spacegroup/spacegroup.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>
#include <cppcrystal/symmetry/primitive.hpp>

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

using namespace cppcrystal;

namespace {

// ---- a corpus spanning the crystal systems ----

Cell primitive_cubic(double a) { // -> Pm-3m
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(Matrix3d::Identity() * a, pos, {0});
}

Cell fcc(double a) { // -> Fm-3m
  Positions pos(4, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  pos.row(2) << 0.5, 0.0, 0.5;
  pos.row(3) << 0.0, 0.5, 0.5;
  return Cell(Matrix3d::Identity() * a, pos, {0, 0, 0, 0});
}

Cell rock_salt(double a) { // -> Fm-3m, two species
  Positions pos(8, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  pos.row(2) << 0.5, 0.0, 0.5;
  pos.row(3) << 0.0, 0.5, 0.5;
  pos.row(4) << 0.5, 0.5, 0.5;
  pos.row(5) << 0.0, 0.0, 0.5;
  pos.row(6) << 0.0, 0.5, 0.0;
  pos.row(7) << 0.5, 0.0, 0.0;
  return Cell(Matrix3d::Identity() * a, pos, {0, 0, 0, 0, 1, 1, 1, 1});
}

Cell bcc_tetragonal() { // -> I4/mmm
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.0;
  l(1, 1) = 4.0;
  l(2, 2) = 6.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(l, pos, {0, 0});
}

Cell rutile() { // -> P4_2/mnm (tetragonal, non-symmorphic)
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
  return Cell(l, pos, {0, 0, 1, 1, 1, 1});
}

Cell orthorhombic_p() { // -> Pmmm
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.0;
  l(1, 1) = 5.0;
  l(2, 2) = 6.0;
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(l, pos, {0});
}

Cell monoclinic_p() { // unique axis b, beta != 90 -> exercises monocli matcher
  Matrix3d l;
  l.col(0) = Vector3d(5.0, 0.0, 0.0);
  l.col(1) = Vector3d(0.0, 6.0, 0.0);
  l.col(2) = Vector3d(1.5, 0.0, 7.0);
  Positions pos(1, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  return Cell(l, pos, {0});
}

Cell hcp() { // -> P6_3/mmc (hexagonal)
  double const a = 3.2;
  double const c = 5.2;
  Matrix3d l;
  l.col(0) = Vector3d(a, 0.0, 0.0);
  l.col(1) = Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, c);
  Positions pos(2, 3);
  pos.row(0) << 1.0 / 3.0, 2.0 / 3.0, 0.25;
  pos.row(1) << 2.0 / 3.0, 1.0 / 3.0, 0.75;
  return Cell(l, pos, {0, 0});
}

Cell rhombohedral_r() { // R-centered hexagonal cell -> trigonal/rhombohedral
  double const a = 3.0;
  double const c = 12.0;
  Matrix3d l;
  l.col(0) = Vector3d(a, 0.0, 0.0);
  l.col(1) = Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, c);
  Positions pos(3, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0;
  pos.row(2) << 2.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0;
  return Cell(l, pos, {0, 0, 0});
}

Cell diamond(double a) { // -> Fd-3m (cubic F, non-symmorphic)
  Positions pos(8, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.0, 0.5, 0.5;
  pos.row(2) << 0.5, 0.0, 0.5;
  pos.row(3) << 0.5, 0.5, 0.0;
  pos.row(4) << 0.25, 0.25, 0.25;
  pos.row(5) << 0.25, 0.75, 0.75;
  pos.row(6) << 0.75, 0.25, 0.75;
  pos.row(7) << 0.75, 0.75, 0.25;
  return Cell(Matrix3d::Identity() * a, pos, {0, 0, 0, 0, 0, 0, 0, 0});
}

Cell zincblende(double a) { // -> F-43m
  Cell d = diamond(a);
  return Cell(d.lattice(), d.positions(), {0, 0, 0, 0, 1, 1, 1, 1});
}

Cell bcc(double a) { // -> Im-3m
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Matrix3d::Identity() * a, pos, {0, 0});
}

Cell cscl(double a) { // -> Pm-3m, two species
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Matrix3d::Identity() * a, pos, {0, 1});
}

Cell wurtzite() { // -> P6_3mc (hexagonal, polar)
  double const a = 3.2;
  double const c = 5.2;
  double const u = 0.375;
  Matrix3d l;
  l.col(0) = Vector3d(a, 0.0, 0.0);
  l.col(1) = Vector3d(-a / 2.0, a * std::sqrt(3.0) / 2.0, 0.0);
  l.col(2) = Vector3d(0.0, 0.0, c);
  Positions pos(4, 3);
  pos.row(0) << 1.0 / 3.0, 2.0 / 3.0, 0.0;
  pos.row(1) << 2.0 / 3.0, 1.0 / 3.0, 0.5;
  pos.row(2) << 1.0 / 3.0, 2.0 / 3.0, u;
  pos.row(3) << 2.0 / 3.0, 1.0 / 3.0, 0.5 + u;
  return Cell(l, pos, {0, 0, 1, 1});
}

Cell c_centered_ortho() { // C-centered orthorhombic -> Cmmm-ish
  Matrix3d l = Matrix3d::Zero();
  l(0, 0) = 4.0;
  l(1, 1) = 6.0;
  l(2, 2) = 7.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  return Cell(l, pos, {0, 0});
}

Cell c_centered_monocli() { // C-centered monoclinic (unique b) -> C2/m
  double const beta = 1.91; // radians, ~109.5 deg
  Matrix3d l;
  l.col(0) = Vector3d(5.0, 0.0, 0.0);
  l.col(1) = Vector3d(0.0, 6.0, 0.0);
  l.col(2) = Vector3d(7.0 * std::cos(beta), 0.0, 7.0 * std::sin(beta));
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.0;
  return Cell(l, pos, {0, 0});
}

Cell triclinic() { // -> P1
  Matrix3d l;
  l.col(0) = Vector3d(4.0, 0.0, 0.0);
  l.col(1) = Vector3d(0.7, 4.3, 0.0);
  l.col(2) = Vector3d(0.5, 0.9, 5.1);
  Positions pos(2, 3);
  pos.row(0) << 0.1, 0.2, 0.3;
  pos.row(1) << 0.6, 0.55, 0.27;
  return Cell(l, pos, {0, 1});
}

struct Case {
  std::string name;
  Cell cell;
};

std::vector<Case> corpus() {
  return {{"primitive_cubic", primitive_cubic(4.0)},
          {"fcc", fcc(4.0)},
          {"rock_salt", rock_salt(5.6)},
          {"bcc_tetragonal", bcc_tetragonal()},
          {"rutile", rutile()},
          {"orthorhombic_p", orthorhombic_p()},
          {"monoclinic_p", monoclinic_p()},
          {"hcp", hcp()},
          {"rhombohedral_r", rhombohedral_r()},
          {"diamond", diamond(3.57)},
          {"zincblende", zincblende(5.43)},
          {"bcc", bcc(3.0)},
          {"cscl", cscl(4.1)},
          {"wurtzite", wurtzite()},
          {"c_centered_ortho", c_centered_ortho()},
          {"c_centered_monocli", c_centered_monocli()},
          {"triclinic", triclinic()}};
}

} // namespace

TEST_CASE("get_dataset number/hall/international match spg_get_dataset",
          "[oracle][spacegroup]") {
  constexpr double symprec = 1e-5;
  for (auto const &[name, cell] : corpus()) {
    INFO("structure: " << name);
    auto const ref = oracle::reference_dataset(cell, symprec);
    REQUIRE(ref.number != 0); // reference must succeed for the comparison

    auto const ours = get_dataset(cell, symprec);
    REQUIRE(ours);

    INFO("ours: SG " << ours->spacegroup_number << " hall "
                     << ours->hall_number << " (" << ours->international_symbol
                     << "); ref: SG " << ref.number << " hall "
                     << ref.hall_number << " (" << ref.international << ")");
    CHECK(ours->spacegroup_number == ref.number);
    CHECK(ours->hall_number == ref.hall_number);
    CHECK(std::string(ours->international_symbol) == ref.international);
    CHECK(std::string(ours->hall_symbol) == ref.hall_symbol);
    CHECK(std::string(ours->choice) == ref.choice);
  }
}

TEST_CASE("standardized lattice / rotation / transformation match the oracle",
          "[oracle][refine]") {
  constexpr double symprec = 1e-5;
  for (auto const &[name, cell] : corpus()) {
    INFO("structure: " << name);
    auto const ref = oracle::reference_dataset(cell, symprec);
    REQUIRE(ref.number != 0);

    auto const prim = symmetry::find_primitive(cell, symprec);
    REQUIRE(prim);
    auto const sg = spacegroup::search_spacegroup(*prim, 0, prim->tolerance,
                                                  prim->angle_tolerance);
    REQUIRE(sg);

    auto const sg2 = refine::find_similar_bravais_lattice(*sg, prim->tolerance);
    Matrix3d const std_lat = refine::conventional_lattice(sg2);
    Matrix3d const std_rot =
        refine::measure_rigid_rotation(sg2.bravais_lattice, std_lat);
    Matrix3d const trans = sg2.bravais_lattice.inverse() * cell.lattice();

    CHECK((std_lat - ref.std_lattice).cwiseAbs().maxCoeff() < 1e-4);
    CHECK((std_rot - ref.std_rotation_matrix).cwiseAbs().maxCoeff() < 1e-4);
    CHECK((trans - ref.transformation_matrix).cwiseAbs().maxCoeff() < 1e-4);
  }
}

namespace {
// True iff every reference std atom has a matching (same type, overlapping)
// atom in ours.
bool same_std_cell(Positions const &ours_pos, std::vector<int> const &ours_types,
                   Positions const &ref_pos, std::vector<int> const &ref_types,
                   Matrix3d const &lattice, double tol) {
  if (ours_types.size() != ref_types.size())
    return false;
  for (Eigen::Index r = 0; r < ref_pos.rows(); ++r) {
    bool matched = false;
    for (Eigen::Index o = 0; o < ours_pos.rows(); ++o) {
      if (ours_types[static_cast<std::size_t>(o)] !=
          ref_types[static_cast<std::size_t>(r)])
        continue;
      if (is_overlap(ours_pos.row(o).transpose(), ref_pos.row(r).transpose(),
                     lattice, tol)) {
        matched = true;
        break;
      }
    }
    if (!matched)
      return false;
  }
  return true;
}
} // namespace

TEST_CASE("get_dataset exposes the full standardized dataset (public API)",
          "[oracle][dataset]") {
  constexpr double symprec = 1e-5;
  for (auto const &[name, cell] : corpus()) {
    INFO("structure: " << name);
    auto const ref = oracle::reference_dataset(cell, symprec);
    REQUIRE(ref.number != 0);

    auto const d = get_dataset(cell, symprec);
    REQUIRE(d);

    // Identity + transformation/standardization matrices.
    CHECK(d->spacegroup_number == ref.number);
    CHECK(d->hall_number == ref.hall_number);
    CHECK((d->transformation_matrix - ref.transformation_matrix)
              .cwiseAbs()
              .maxCoeff() < 1e-4);
    CHECK((d->std_lattice - ref.std_lattice).cwiseAbs().maxCoeff() < 1e-4);
    CHECK((d->std_rotation_matrix - ref.std_rotation_matrix)
              .cwiseAbs()
              .maxCoeff() < 1e-4);

    // Standardized cell + per-atom Wyckoff data.
    CHECK(static_cast<int>(d->std_positions.rows()) == ref.n_std_atoms);
    CHECK(same_std_cell(d->std_positions, d->std_types, ref.std_positions,
                        ref.std_types, d->std_lattice, 1e-4));
    CHECK(d->wyckoffs == ref.wyckoffs);
    CHECK(d->equivalent_atoms == ref.equivalent_atoms);
    CHECK(d->site_symmetry_symbols.size() == ref.site_symmetry_symbols.size());
    for (std::size_t i = 0; i < d->site_symmetry_symbols.size(); ++i) {
      INFO("atom " << i);
      CHECK(d->site_symmetry_symbols[i] == ref.site_symmetry_symbols[i]);
    }
  }
}

TEST_CASE("standardized cell + Wyckoffs + equivalent atoms match the oracle",
          "[oracle][wyckoff]") {
  constexpr double symprec = 1e-5;
  for (auto const &[name, cell] : corpus()) {
    INFO("structure: " << name);
    auto const ref = oracle::reference_dataset(cell, symprec);
    REQUIRE(ref.number != 0);

    auto const prim = symmetry::find_primitive(cell, symprec);
    REQUIRE(prim);
    auto const sg = spacegroup::search_spacegroup(*prim, 0, prim->tolerance,
                                                  prim->angle_tolerance);
    REQUIRE(sg);
    auto const ops = symmetry::find_symmetry(cell, prim->tolerance);
    REQUIRE(ops);

    auto const sg2 = refine::find_similar_bravais_lattice(*sg, prim->tolerance);
    auto const std = refine::get_wyckoff_positions(
        sg2, prim->cell, cell, *ops, prim->mapping_table, prim->tolerance);
    REQUIRE(std);

    // Standardized cell: same atom count, same atoms (set, mod lattice).
    CHECK(static_cast<int>(std->bravais.size()) == ref.n_std_atoms);
    CHECK(same_std_cell(std->bravais.positions(), std->bravais.types(),
                        ref.std_positions, ref.std_types,
                        std->bravais.lattice(), 1e-4));

    // Per-atom Wyckoff letters, site-symmetry symbols, equivalent atoms.
    CHECK(std->wyckoffs == ref.wyckoffs);
    CHECK(std->equivalent_atoms == ref.equivalent_atoms);
    CHECK(std->site_symmetry_symbols.size() ==
          ref.site_symmetry_symbols.size());
    for (std::size_t i = 0; i < std->site_symmetry_symbols.size(); ++i) {
      INFO("atom " << i);
      CHECK(std->site_symmetry_symbols[i] == ref.site_symmetry_symbols[i]);
    }
  }
}

TEST_CASE("get_dataset operations equal the input cell symmetry set",
          "[oracle][spacegroup]") {
  constexpr double symprec = 1e-5;
  for (auto const &[name, cell] : corpus()) {
    INFO("structure: " << name);
    auto const ref = oracle::reference_dataset(cell, symprec);
    REQUIRE(ref.number != 0);

    auto const ours = get_dataset(cell, symprec);
    REQUIRE(ours);
    CHECK(static_cast<int>(ours->operations.size()) == ref.n_operations);

    // Set-equal to the reference operations (rotations exact, translations mod
    // lattice). spglib orders them differently than we discover them.
    auto same_set = [](SymmetryOperations const &a, SymmetryOperations const &b) {
      if (a.size() != b.size())
        return false;
      for (auto const &oa : a) {
        bool matched = false;
        for (auto const &ob : b)
          if (same_operation(oa, ob, 1e-4)) {
            matched = true;
            break;
          }
        if (!matched)
          return false;
      }
      return true;
    };
    CHECK(same_set(ours->operations, ref.operations));
  }
}
