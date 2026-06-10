#include <spglib/dataset.hpp>

#include <spglib/refine/operations.hpp>
#include <spglib/refine/refinement.hpp>
#include <spglib/refine/standardize.hpp>
#include <spglib/spacegroup/spacegroup.hpp>
#include <spglib/symmetry/pointgroup.hpp>
#include <spglib/symmetry/primitive.hpp>

// Port of the determination + refinement pipeline (determination.c
// det_determine_all + spglib.c set_dataset), 3D space-group path. The outer
// loop progressively tightens the tolerance (REDUCE_RATE_OUTER) in case a given
// tolerance yields an inconsistent cell. Operations are taken from the input
// cell's symmetry search (set-equal to spglib's exact refined operations); the
// idealized refined operation values are a later refinement.
namespace spglib {

namespace {

constexpr int kNumAttemptOuter = 10;     // determination.c NUM_ATTEMPT_OUTER
constexpr double kReduceRateOuter = 0.9; // determination.c REDUCE_RATE_OUTER

} // namespace

Result<Dataset> get_dataset(Cell const &cell, double symprec,
                            AngleTolerance angle_tolerance, int hall_number) {
  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttemptOuter;
       ++attempt, tolerance *= kReduceRateOuter) {
    auto const primitive =
        symmetry::find_primitive(cell, tolerance, angle_tolerance);
    if (!primitive) {
      continue;
    }
    double const tol = primitive->tolerance;
    auto const sg = spacegroup::search_spacegroup(*primitive, hall_number, tol,
                                                  primitive->angle_tolerance);
    if (!sg) {
      continue;
    }

    // Standardize: orient the bravais lattice, recover the exact operations in
    // the input cell, and build the idealized cell + Wyckoff data.
    spacegroup::Spacegroup const sg2 =
        refine::find_similar_bravais_lattice(*sg, tol);

    auto const operations =
        refine::refined_operations(sg2, primitive->cell, cell, tol);
    if (!operations) {
      continue;
    }

    auto const std = refine::get_wyckoff_positions(
        sg2, primitive->cell, cell, *operations, primitive->mapping_table, tol);
    if (!std) {
      continue;
    }

    data::SpacegroupType const &t = sg2.type;
    Matrix3d const std_lattice = std->bravais.lattice();
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
        .transformation_matrix = sg2.bravais_lattice.inverse() * cell.lattice(),
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
            refine::measure_rigid_rotation(sg2.bravais_lattice, std_lattice),
        .std_mapping_to_primitive = std->std_mapping_to_primitive,

        .primitive_lattice = primitive->cell.lattice(),
        .mapping_to_primitive = primitive->mapping_table,
    };
  }
  return leaf::new_error(e_spacegroup_search_failed{});
}

} // namespace spglib
