#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>

#include <utility>

namespace cppcrystal::analysis {

MagneticSymmetryAnalyzer
MagneticSymmetryAnalyzer::from_cell(MagneticCell cell, MagneticTolerance tol) {
  return MagneticSymmetryAnalyzer{std::move(cell), tol};
}

Result<MagneticDataset const *>
MagneticSymmetryAnalyzer::cached_dataset() const {
  return dataset_.get([&] { return get_magnetic_dataset(cell_, tol_); });
}

Result<MagneticCell> MagneticSymmetryAnalyzer::standardized_cell() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return MagneticCell{
      Cell{Lattice{ds->std_lattice}, ds->std_positions, ds->std_types},
      ds->std_tensors, cell_.kind()};
}

Result<void> MagneticSymmetryAnalyzer::warm() const {
  BOOST_LEAF_CHECK(cached_dataset());
  return {};
}

} // namespace cppcrystal::analysis
