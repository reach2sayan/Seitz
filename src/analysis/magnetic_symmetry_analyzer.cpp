#include <spglib/analysis/magnetic_symmetry_analyzer.hpp>

#include <spglib/symmetry/find_symmetry.hpp>

#include <utility>

namespace spglib::analysis {

MagneticSymmetryAnalyzer MagneticSymmetryAnalyzer::from_cell(
    MagneticCell cell, bool is_axial, double symprec,
    AngleTolerance angle_tolerance, std::optional<double> mag_symprec) {
  return MagneticSymmetryAnalyzer{std::move(cell), is_axial,
                                  Tolerance{symprec, angle_tolerance},
                                  mag_symprec};
}

Result<MagneticDataset const *> MagneticSymmetryAnalyzer::cached_dataset()
    const {
  if (!dataset_) {
    BOOST_LEAF_AUTO(ds, get_magnetic_dataset(cell_, is_axial_, tol_.symprec,
                                             tol_.angle_tolerance,
                                             mag_symprec_));
    dataset_ = std::move(ds);
  }
  return &*dataset_;
}

Result<MagneticDataset> MagneticSymmetryAnalyzer::dataset() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return *ds;
}

Result<MagneticSymmetryOperations> MagneticSymmetryAnalyzer::operations()
    const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->operations;
}

Result<spin::MagneticSymmetrySearch>
MagneticSymmetryAnalyzer::symmetry_search() const {
  if (!symmetry_search_) {
    BOOST_LEAF_AUTO(sym_nonspin,
                    symmetry::find_symmetry(cell_.cell(), tol_.symprec,
                                            tol_.angle_tolerance));
    BOOST_LEAF_AUTO(search, spin::operations_with_site_tensors(
                                sym_nonspin, cell_, /*with_time_reversal=*/true,
                                is_axial_, tol_.symprec, tol_.angle_tolerance,
                                mag_symprec_));
    symmetry_search_ = std::move(search);
  }
  return *symmetry_search_;
}

Result<int> MagneticSymmetryAnalyzer::uni_number() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->uni_number;
}

Result<int> MagneticSymmetryAnalyzer::hall_number() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->hall_number;
}

Result<data::MagneticSpacegroupType>
MagneticSymmetryAnalyzer::spacegroup_type() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return data::magnetic_spacegroup_type(ds->uni_number);
}

Result<std::vector<int>> MagneticSymmetryAnalyzer::equivalent_atoms() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->equivalent_atoms;
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

} // namespace spglib::analysis
