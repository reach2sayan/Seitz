#include <cppcrystal/dataset.hpp>

#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/validation.hpp>
#include <cppcrystal/refine/operations.hpp>
#include <cppcrystal/refine/refinement.hpp>
#include <cppcrystal/refine/standardize.hpp>
#include <cppcrystal/spacegroup/spacegroup.hpp>
#include <cppcrystal/symmetry/pointgroup.hpp>
#include <cppcrystal/symmetry/primitive.hpp>

// The determination + refinement pipeline, 3D space-group path. The outer loop
// progressively tightens the tolerance in case a given tolerance yields an
// inconsistent cell. Operations are taken from the input cell's symmetry search;
// the idealized refined operation values are a later refinement.
namespace cppcrystal {

namespace {

constexpr int kNumAttemptOuter = 10;
constexpr double kReduceRateOuter = 0.9; // tolerance shrink factor per attempt

} // namespace

Result<Dataset> get_dataset(Cell const &cell, Tolerance const &tol,
                            int hall_number) {
  if (auto valid = validate_cell(cell); !valid) {
    return valid.error();
  }
  Tolerance attempt_tol = tol;
  for (int attempt = 0; attempt < kNumAttemptOuter;
       ++attempt, attempt_tol.symprec *= kReduceRateOuter) {
    auto const primitive = symmetry::find_primitive(cell, attempt_tol);
    if (!primitive) {
      continue;
    }
    Tolerance const &found = primitive->tolerance;
    auto const sg =
        spacegroup::search_spacegroup(*primitive, hall_number, found);
    if (!sg) {
      continue;
    }

    // Standardize: orient the bravais lattice, recover the exact operations in
    // the input cell, and build the idealized cell + Wyckoff data.
    spacegroup::Spacegroup const sg2 =
        refine::find_similar_bravais_lattice(*sg, found.symprec);

    auto const operations =
        refine::refined_operations(sg2, primitive->cell, cell, found.symprec);
    if (!operations) {
      continue;
    }

    auto const std =
        refine::get_wyckoff_positions(sg2, primitive->cell, cell, *operations,
                                      primitive->mapping_table, found.symprec);
    if (!std) {
      continue;
    }

    data::SpacegroupType const &t = sg2.type;
    Matrix3d const std_lattice = std->bravais.lattice().matrix();
    return Dataset{
        .spacegroup_number = t.number,
        .hall_number = t.hall_number,
        .international_symbol = t.international_short,
        .hall_symbol = t.hall_symbol,
        .choice = t.choice,
        .pointgroup_number = t.pointgroup_number,
        .pointgroup_symbol =
            symmetry::pointgroup_by_number(t.pointgroup_number).symbol,

        .bravais_lattice = sg2.bravais_lattice,
        .transformation_matrix = sg2.bravais_lattice.inverse() * cell.lattice().matrix(),
        .origin_shift = sg2.origin_shift,
        .operations = *operations,

        .wyckoffs = std->wyckoffs,
        .site_symmetry_symbols = std->site_symmetry_symbols,
        .equivalent_atoms = std->equivalent_atoms,
        .crystallographic_orbits = std->crystallographic_orbits,

        .std_lattice = std_lattice,
        .std_positions = std->bravais.positions(),
        .std_types = std->bravais.types(),
        .std_rotation_matrix =
            Lattice{sg2.bravais_lattice}.rigid_rotation_to(Lattice{std_lattice}),
        .std_mapping_to_primitive = std->std_mapping_to_primitive,

        .primitive_lattice = primitive->cell.lattice().matrix(),
        .mapping_to_primitive = primitive->mapping_table,
        .aperiodic_axis = aperiodic_axis(cell.periodicity()),
    };
  }
  return leaf::new_error(e_spacegroup_search_failed{});
}

} // namespace cppcrystal
