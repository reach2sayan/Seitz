#include <cppcrystal/dataset.hpp>

#include <cppcrystal/core/lattice.hpp>
#include "core/validation.hpp"
#include "refine/refinement.hpp"
#include "refine/standardize.hpp"
#include "spacegroup/spacegroup.hpp"
#include "symmetry/pointgroup.hpp"
#include "core/family_dispatch.hpp"
#include "symmetry/primitive.hpp"

// The determination + refinement pipeline, 3D space-group path. The outer loop
// progressively tightens the tolerance in case a given tolerance yields an
// inconsistent cell. Operations are taken from the input cell's symmetry search;
// the idealized refined operation values are a later refinement.
namespace cppcrystal {

namespace {

constexpr int kNumAttemptOuter = 10;
constexpr double kReduceRateOuter = 0.9; // tolerance shrink factor per attempt

// The determination + refinement pipeline for one family. Everything below the
// dispatch in get_dataset is compile-time-specialised on it.
template <GroupFamily F>
[[nodiscard]] Result<Dataset> run_pipeline(Cell const &cell,
                                           Tolerance const &tol,
                                           std::optional<HallNumber> setting) {
  Tolerance attempt_tol = tol;
  for (int attempt = 0; attempt < kNumAttemptOuter;
       ++attempt, attempt_tol.symprec *= kReduceRateOuter) {
    symmetry::PrimitiveFinder<F> const finder(cell, attempt_tol);
    auto const primitive = finder.find();
    if (!primitive) {
      continue;
    }
    Tolerance const &found = primitive->tolerance;
    spacegroup::SpacegroupMatcher<F> const matcher(*primitive, setting);
    auto const sg = matcher.search();
    if (!sg) {
      continue;
    }

    // Standardize: orient the bravais lattice, recover the exact operations in
    // the input cell, and build the idealized cell + Wyckoff data.
    refine::Refinement<F> const refinement =
        refine::Refinement<F>{*sg, primitive->cell, cell, found}
            .similar_bravais();
    SpacegroupMatch const &sg2 = refinement.matched();

    auto const operations = refinement.operations();
    if (!operations) {
      continue;
    }

    auto const std =
        refinement.standardize(*operations, primitive->mapping_table);
    if (!std) {
      continue;
    }

    data::SpacegroupType const &t = sg2.type();
    Matrix3d const std_lattice = std->bravais.lattice().matrix();
    return Dataset{
        .spacegroup_number = t.number,
        .hall = sg2.hall,
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

} // namespace

Result<Dataset> get_dataset(Cell const &cell, Tolerance const &tol,
                            std::optional<HallNumber> setting) {
  if (auto valid = validate_cell(cell); !valid) {
    return valid.error();
  }
  // The one runtime branch on the family; the pipeline below it is templated.
  return dispatch_family(cell.periodicity(), [&]<GroupFamily F>() {
    return run_pipeline<F>(cell, tol, setting);
  });
}

} // namespace cppcrystal
