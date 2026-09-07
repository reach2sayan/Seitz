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

  // 2. Identify the magnetic space-group type: the public door.
  BOOST_LEAF_AUTO(match,
                  search.operations.spacegroup(cell_.cell().lattice(), tol_));

  // 3. Idealize positions and site tensors, then transform into the
  //    standardized setting.
  MagneticCell const exact = spin_search.idealized<kTimeReversal>(search);
  magnetic::MagneticIdentification const identification(
      cell_.cell().lattice(), search.operations, tol_);
  BOOST_LEAF_AUTO(standardized, identification.transform(exact, match));

  // Base first, then the cell-level products. `search` is dead after this
  // aggregate except for primitive_lattice, which is a different member.
  return MagneticDataset{std::move(match), std::move(search.operations),
                         std::move(search.equivalent_atoms),
                         std::move(standardized),
                         Lattice{search.primitive_lattice}};
}

} // namespace seitz::analysis
