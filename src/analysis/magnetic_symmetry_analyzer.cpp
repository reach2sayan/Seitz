#include <seitz/analysis/magnetic_symmetry_analyzer.hpp>

#include "core/validation.hpp"
#include "magnetic/identify.hpp"
#include "spin/search.hpp"
#include "symmetry/search.hpp"

#include <utility>

namespace seitz::analysis {

Result<MagneticDataset> MagneticSymmetryAnalyzer::determine() const {
  // The magnetic determination always searches with time reversal (the full
  // family space group), and is a 3D path only.
  constexpr TimeReversal kTimeReversal = TimeReversal::on;
  if (auto valid = validate_cell(cell_.cell()); !valid) {
    return valid.error();
  }

  // 1. Magnetic symmetry of the input cell.
  symmetry::SymmetrySearch<GroupFamily::space> const spatial(cell_.cell(),
                                                             tol_);
  BOOST_LEAF_AUTO(sym_nonspin, spatial.operations());
  spin::SpinSearch const spin_search(cell_, sym_nonspin, tol_);
  BOOST_LEAF_AUTO(search, spin_search.operations<kTimeReversal>());

  // 2. Identify the magnetic space-group type.
  magnetic::MagneticIdentification const identification(
      cell_.cell().lattice(), search.operations, tol_);
  BOOST_LEAF_AUTO(ident, identification.identify());

  // 3. Idealize positions and site tensors, then transform into the
  //    standardized setting.
  MagneticCell const exact = spin_search.idealized<kTimeReversal>(search);
  BOOST_LEAF_AUTO(standardized, identification.transform(exact, ident));

  return MagneticDataset{
      .uni = ident.uni,
      .type = ident.msg_type,
      .hall = ident.hall,
      .setting = {.transformation = ident.transformation_matrix,
                  .origin_shift = ident.origin_shift,
                  .rigid_rotation = ident.std_rotation_matrix},
      // Moved: `search` is dead after this aggregate except for
      // primitive_lattice below, which is a different member.
      .operations = std::move(search.operations),
      .equivalent_atoms = std::move(search.equivalent_atoms),
      .standardized = std::move(standardized),
      .primitive = Lattice{search.primitive_lattice},
  };
}

} // namespace seitz::analysis
