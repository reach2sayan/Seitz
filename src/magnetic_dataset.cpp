#include <cppcrystal/magnetic_dataset.hpp>

#include <cppcrystal/core/validation.hpp>
#include <cppcrystal/magnetic/magnetic_spacegroup.hpp>
#include <cppcrystal/spin/spin.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>

namespace cppcrystal {

Result<MagneticDataset> get_magnetic_dataset(MagneticCell const &cell,
                                             MagneticTolerance const &tol) {
  // The magnetic dataset always searches with time reversal (the full family
  // space group).
  constexpr TimeReversal kTimeReversal = TimeReversal::on;
  if (auto valid = validate_cell(cell.cell()); !valid) {
    return valid.error();
  }
  double const symprec = tol.symprec;

  // 1. Magnetic symmetry of the input cell.
  BOOST_LEAF_AUTO(sym_nonspin, symmetry::find_symmetry(cell.cell(), tol));
  BOOST_LEAF_AUTO(search, spin::operations_with_site_tensors(
                              sym_nonspin, cell, kTimeReversal, tol));

  // 2. Identify the magnetic space-group type.
  BOOST_LEAF_AUTO(ident, magnetic::identify_magnetic_spacegroup_type(
                             cell.cell().lattice().matrix(), search.operations,
                             symprec));

  // 3. Idealize positions and site tensors.
  MagneticCell const exact = spin::idealized_cell(search, cell, kTimeReversal);

  // 4. Transform the idealized cell into the standardized setting.
  BOOST_LEAF_AUTO(std_cell,
                  magnetic::transform_cell(
                      exact, ident.transformation_matrix, ident.origin_shift,
                      ident.std_rotation_matrix, search.operations, symprec));

  MagneticDataset const dataset{
      .uni_number = ident.uni_number,
      .msg_type = ident.msg_type,
      .hall_number = ident.hall_number,
      .tensor_rank = cell.rank(),

      .operations = search.operations,
      .equivalent_atoms = search.equivalent_atoms,

      .transformation_matrix = ident.transformation_matrix,
      .origin_shift = ident.origin_shift,

      .std_lattice = std_cell.cell().lattice().matrix(),
      .std_positions = std_cell.cell().positions(),
      .std_types = std_cell.cell().types(),
      .std_tensors = std_cell.tensors(),
      .std_rotation_matrix = ident.std_rotation_matrix,

      .primitive_lattice = search.primitive_lattice,
  };
  return dataset;
}

} // namespace cppcrystal
