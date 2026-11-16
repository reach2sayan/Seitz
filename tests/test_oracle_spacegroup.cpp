#include "oracle.hpp"
#include <seitz/core/operation_set.hpp>

#include "core/overlap.hpp"
#include "refine/refinement.hpp"
#include "refine/standardize.hpp"
#include "spacegroup/spacegroup.hpp"
#include "symmetry/primitive.hpp"
#include "symmetry/search.hpp"
#include <seitz/analysis/symmetry_analyzer.hpp>

#include "helpers.hpp"
#include "structures.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <string>
#include <utility>
#include <vector>

using namespace seitz;
using seitz::test::space_hall;

namespace {

using namespace seitz::test;

constexpr double kSymprec = 1e-5;

// One structure with its reference dataset, computed once for every case.
struct Case {
  std::string name;
  Cell cell;
  oracle::RefDataset ref;
};

std::vector<Case> const &corpus() {
  static std::vector<Case> const cases = [] {
    std::vector<std::pair<std::string, Cell>> const cells{
        {"primitive_cubic", primitive_cubic(4.0)},
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
        {"cscl", rocksalt_motif(4.1)},
        {"wurtzite", wurtzite()},
        {"c_centered_ortho", c_centered_ortho()},
        {"c_centered_monocli", c_centered_monocli()},
        {"triclinic", triclinic()}};
    return std::vector<Case>{
        std::from_range, cells | std::views::transform([](auto const &c) {
                           return Case{c.first, c.second,
                                       oracle::reference_dataset(c.second,
                                                                 kSymprec)};
                         })};
  }();
  return cases;
}

// The pipeline stage by stage, up to the space-group match the refinement
// starts from.
struct Matched {
  symmetry::Primitive prim;
  SpacegroupMatch sg;
};
Matched matched(Cell const &cell) {
  auto const prim =
      symmetry::PrimitiveFinder<GroupFamily::space>{cell, {kSymprec}}.find();
  REQUIRE(prim);
  spacegroup::SpacegroupMatcher<GroupFamily::space> const matcher(
      *prim, std::nullopt);
  auto const sg = matcher.search();
  REQUIRE(sg);
  return {*prim,
          refine::find_similar_bravais_lattice(*sg, prim->tolerance.symprec)};
}

// True iff every reference std atom has a matching (same type, overlapping)
// atom in ours.
bool same_std_cell(Cell const &ours, oracle::RefDataset const &ref,
                   double tol) {
  if (ours.types().size() != ref.std_types.size()) {
    return false;
  }
  return std::ranges::all_of(
      std::views::enumerate(ref.std_types), [&](auto const &entry) {
        auto const [r, ref_type] = entry;
        Vector3d const ref_pos = ref.std_positions.row(r).transpose();
        return std::ranges::any_of(ours.atoms(), [&](auto const &atom) {
          auto const &[pos, type] = atom;
          return type == ref_type &&
                 coincident(pos, ref_pos, ours.lattice().matrix(), tol,
                            all_periodic());
        });
      });
}

} // namespace

TEST_CASE("get_dataset number/hall/international match spg_get_dataset",
          "[oracle][spacegroup]") {
  for (auto const &[name, cell, ref] : corpus()) {
    INFO("structure: " << name);
    REQUIRE(ref.number != 0); // reference must succeed for the comparison

    auto const ours = test::dataset_of(cell, {kSymprec});
    REQUIRE(ours);

    INFO("ours: SG " << data::spacegroup_type(ours->hall).number << " hall "
                     << ours->hall.index() << " ("
                     << data::spacegroup_type(ours->hall).international_short
                     << "); ref: SG " << ref.number << " hall "
                     << ref.hall_number << " (" << ref.international << ")");
    CHECK(data::spacegroup_type(ours->hall).number == ref.number);
    CHECK(ours->hall == space_hall(ref.hall_number));
    CHECK(std::string(data::spacegroup_type(ours->hall).international_short) ==
          ref.international);
    CHECK(std::string(data::spacegroup_type(ours->hall).hall_symbol) ==
          ref.hall_symbol);
    CHECK(std::string(data::spacegroup_type(ours->hall).choice) == ref.choice);
  }
}

TEST_CASE("standardized lattice / rotation / transformation match the oracle",
          "[oracle][refine]") {
  for (auto const &[name, cell, ref] : corpus()) {
    INFO("structure: " << name);
    REQUIRE(ref.number != 0);

    auto const [prim, sg2] = matched(cell);
    Lattice const std_lattice = refine::conventional_lattice(sg2);
    Matrix3d const std_lat = std_lattice.matrix();
    Matrix3d const std_rot =
        Lattice{sg2.bravais_lattice}.rigid_rotation_to(std_lattice);
    Matrix3d const trans =
        sg2.bravais_lattice.inverse() * cell.lattice().matrix();

    CHECK((std_lat - ref.std_lattice).cwiseAbs().maxCoeff() < 1e-4);
    CHECK((std_rot - ref.std_rotation_matrix).cwiseAbs().maxCoeff() < 1e-4);
    CHECK((trans - ref.transformation_matrix).cwiseAbs().maxCoeff() < 1e-4);
  }
}

TEST_CASE("get_dataset exposes the full standardized dataset (public API)",
          "[oracle][dataset]") {
  for (auto const &[name, cell, ref] : corpus()) {
    INFO("structure: " << name);
    REQUIRE(ref.number != 0);

    auto const d = test::dataset_of(cell, {kSymprec});
    REQUIRE(d);

    // Identity + transformation/standardization matrices.
    CHECK(data::spacegroup_type(d->hall).number == ref.number);
    CHECK(d->hall == space_hall(ref.hall_number));
    CHECK((d->setting.transformation - ref.transformation_matrix)
              .cwiseAbs()
              .maxCoeff() < 1e-4);
    CHECK((d->standardized.lattice().matrix() - ref.std_lattice)
              .cwiseAbs()
              .maxCoeff() < 1e-4);
    CHECK((d->setting.rigid_rotation - ref.std_rotation_matrix)
              .cwiseAbs()
              .maxCoeff() < 1e-4);

    // Standardized cell + per-atom Wyckoff data.
    CHECK(static_cast<int>(d->standardized.positions().rows()) ==
          ref.n_std_atoms);
    CHECK(same_std_cell(d->standardized, ref, 1e-4));
    oracle::check_sites(d->sites, ref);
  }
}

TEST_CASE("standardized cell + Wyckoffs + equivalent atoms match the oracle",
          "[oracle][wyckoff]") {
  for (auto const &[name, cell, ref] : corpus()) {
    INFO("structure: " << name);
    REQUIRE(ref.number != 0);

    auto const [prim, sg2] = matched(cell);
    auto const ops =
        symmetry::SymmetrySearch<GroupFamily::space>{cell, prim.tolerance}
            .operations();
    REQUIRE(ops);

    refine::Refinement<GroupFamily::space> const refinement(
        sg2, prim.cell, cell, prim.tolerance);
    auto const std = refinement.standardize(*ops, prim.mapping_table);
    REQUIRE(std);

    // Standardized cell: same atom count, same atoms (set, mod lattice).
    CHECK(static_cast<int>(std->bravais.size()) == ref.n_std_atoms);
    CHECK(same_std_cell(std->bravais, ref, 1e-4));

    // Per-atom Wyckoff letters, site-symmetry symbols, equivalent atoms.
    oracle::check_sites(std->sites, ref);
  }
}

TEST_CASE("get_dataset operations equal the input cell symmetry set",
          "[oracle][spacegroup]") {
  for (auto const &[name, cell, ref] : corpus()) {
    INFO("structure: " << name);
    REQUIRE(ref.number != 0);

    auto const ours = test::dataset_of(cell, {kSymprec});
    REQUIRE(ours);
    CHECK(static_cast<int>(ours->operations.size()) == ref.n_operations);

    // Set-equal to the reference operations (rotations exact, translations mod
    // lattice). spglib orders them differently than we discover them.
    CHECK(same_operation_set(ours->operations, ref.operations, 1e-4));
  }
}
