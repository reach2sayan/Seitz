#include <cppcrystal/magnetic_dataset.hpp>

#include <cppcrystal/core/validation.hpp>
#include <cppcrystal/magnetic/magnetic_spacegroup.hpp>
#include <cppcrystal/spin/spin.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>

namespace cppcrystal {

Result<MagneticDataset>
get_magnetic_dataset(MagneticCell const &cell, bool is_axial, double symprec,
                     AngleTolerance angle_tolerance,
                     std::optional<double> mag_symprec) {
  // The magnetic dataset always searches with time reversal (the full family
  // space group).
  constexpr bool kWithTimeReversal = true;
  if (auto valid = validate_cell(cell.cell()); !valid) {
    return valid.error();
  }

  // 1. Magnetic symmetry of the input cell.
  BOOST_LEAF_AUTO(sym_nonspin, symmetry::find_symmetry(cell.cell(), symprec,
                                                       angle_tolerance));
  BOOST_LEAF_AUTO(search, spin::operations_with_site_tensors(
                              sym_nonspin, cell, kWithTimeReversal, is_axial,
                              symprec, angle_tolerance, mag_symprec));

  // 2. Identify the magnetic space-group type.
  BOOST_LEAF_AUTO(ident,
                  magnetic::identify_magnetic_spacegroup_type(
                      cell.cell().lattice(), search.operations, symprec));

  // 3. Idealize positions and site tensors.
  MagneticCell const exact =
      spin::idealized_cell(search, cell, kWithTimeReversal, is_axial);

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

      .std_lattice = std_cell.cell().lattice(),
      .std_positions = std_cell.cell().positions(),
      .std_types = std_cell.cell().types(),
      .std_tensors = std_cell.tensors(),
      .std_rotation_matrix = ident.std_rotation_matrix,

      .primitive_lattice = search.primitive_lattice,
  };
  return dataset;
}

} // namespace cppcrystal
