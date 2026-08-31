#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>

#include <cppcrystal/symmetry/find_symmetry.hpp>

#include <utility>

namespace cppcrystal::analysis {

MagneticSymmetryAnalyzer MagneticSymmetryAnalyzer::from_cell(
    MagneticCell cell, bool is_axial, double symprec,
    AngleTolerance angle_tolerance, std::optional<double> mag_symprec) {
  return MagneticSymmetryAnalyzer{std::move(cell), is_axial,
                                  Tolerance{symprec, angle_tolerance},
                                  mag_symprec};
}

Result<MagneticDataset const *>
MagneticSymmetryAnalyzer::cached_dataset() const {
  return dataset_.get([&] {
    return get_magnetic_dataset(cell_, is_axial_, tol_.symprec,
                                tol_.angle_tolerance, mag_symprec_);
  });
}

Result<spin::MagneticSymmetrySearch>
MagneticSymmetryAnalyzer::symmetry_search() const {
  BOOST_LEAF_AUTO(
      search,
      symmetry_search_.get([&]() -> Result<spin::MagneticSymmetrySearch> {
        BOOST_LEAF_AUTO(sym_nonspin,
                        symmetry::find_symmetry(cell_.cell(), tol_.symprec,
                                                tol_.angle_tolerance));
        return spin::operations_with_site_tensors(
            sym_nonspin, cell_, /*with_time_reversal=*/true, is_axial_,
            tol_.symprec, tol_.angle_tolerance, mag_symprec_);
      }));
  return *search;
}

Result<MagneticCell> MagneticSymmetryAnalyzer::standardized_cell() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return MagneticCell{Cell{ds->std_lattice, ds->std_positions, ds->std_types},
                      ds->std_tensors};
}

Result<void> MagneticSymmetryAnalyzer::warm() const {
  BOOST_LEAF_CHECK(cached_dataset());
  BOOST_LEAF_CHECK(symmetry_search());
  return {};
}

} // namespace cppcrystal::analysis
