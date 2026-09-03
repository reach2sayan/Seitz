#include <cppcrystal/magnetic_dataset.hpp>

#include "core/validation.hpp"
#include "magnetic/identify.hpp"
#include "spin/search.hpp"
#include "symmetry/search.hpp"

namespace cppcrystal {

Result<MagneticDataset> get_magnetic_dataset(MagneticCell const &cell,
                                             MagneticTolerance const &tol) {
  // The magnetic dataset always searches with time reversal (the full family
  // space group).
  constexpr TimeReversal kTimeReversal = TimeReversal::on;
  if (auto valid = validate_cell(cell.cell()); !valid) {
    return valid.error();
  }

  // 1. Magnetic symmetry of the input cell.
  // The magnetic search is a 3D path only.
  symmetry::SymmetrySearch<GroupFamily::space> const spatial(cell.cell(), tol);
  BOOST_LEAF_AUTO(sym_nonspin, spatial.operations());
  spin::SpinSearch const spin_search(cell, sym_nonspin, tol);
  BOOST_LEAF_AUTO(search, spin_search.operations<kTimeReversal>());

  // 2. Identify the magnetic space-group type.
  magnetic::MagneticIdentification const identification(
      cell.cell().lattice(), search.operations, tol);
  BOOST_LEAF_AUTO(ident, identification.identify());

  // 3. Idealize positions and site tensors.
  MagneticCell const exact = spin_search.idealized<kTimeReversal>(search);

  // 4. Transform the idealized cell into the standardized setting.
  BOOST_LEAF_AUTO(std_cell, identification.transform(exact, ident));

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
