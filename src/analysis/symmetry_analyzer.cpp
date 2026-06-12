#include <spglib/analysis/symmetry_analyzer.hpp>

#include <utility>

namespace spglib::analysis {

SymmetryAnalyzer SymmetryAnalyzer::from_cell(Cell cell, double symprec,
                                             AngleTolerance angle_tolerance,
                                             int hall_number) {
  return SymmetryAnalyzer{std::move(cell), Tolerance{symprec, angle_tolerance},
                          hall_number};
}

SymmetryAnalyzer SymmetryAnalyzer::from_layer_cell(
    Cell cell, int aperiodic_axis, double symprec,
    AngleTolerance angle_tolerance) {
  cell.set_aperiodic_axis(aperiodic_axis);
  return SymmetryAnalyzer{std::move(cell), Tolerance{symprec, angle_tolerance},
                          /*hall_number=*/0};
}

Result<symmetry::Primitive const *> SymmetryAnalyzer::cached_primitive() const {
  if (!primitive_) {
    BOOST_LEAF_AUTO(prim, symmetry::find_primitive(cell_, tol_.symprec,
                                                   tol_.angle_tolerance));
    primitive_ = std::move(prim);
  }
  return &*primitive_;
}

Result<spacegroup::Spacegroup const *>
SymmetryAnalyzer::cached_spacegroup() const {
  if (!spacegroup_) {
    BOOST_LEAF_AUTO(prim, cached_primitive());
    BOOST_LEAF_AUTO(sg, spacegroup::search_spacegroup(
                            *prim, hall_number_, tol_.symprec,
                            tol_.angle_tolerance));
    spacegroup_ = std::move(sg);
  }
  return &*spacegroup_;
}

Result<Dataset const *> SymmetryAnalyzer::cached_dataset() const {
  if (!dataset_) {
    BOOST_LEAF_AUTO(ds, get_dataset(cell_, tol_.symprec, tol_.angle_tolerance,
                                    hall_number_));
    dataset_ = std::move(ds);
  }
  return &*dataset_;
}

Result<Dataset> SymmetryAnalyzer::dataset() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return *ds;
}

Result<SymmetryOperations> SymmetryAnalyzer::operations() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->operations;
}

Result<SymmetryOperations> SymmetryAnalyzer::cell_operations() const {
  if (!cell_operations_) {
    BOOST_LEAF_AUTO(ops, symmetry::find_symmetry(cell_, tol_.symprec,
                                                 tol_.angle_tolerance));
    cell_operations_ = std::move(ops);
  }
  return *cell_operations_;
}

Result<PointSymmetry> SymmetryAnalyzer::lattice_symmetry() const {
  if (!lattice_symmetry_) {
    BOOST_LEAF_AUTO(ps, symmetry::lattice_symmetry(cell_, tol_.symprec,
                                                   tol_.angle_tolerance));
    lattice_symmetry_ = std::move(ps);
  }
  return *lattice_symmetry_;
}

Result<data::SpacegroupType> SymmetryAnalyzer::spacegroup_type() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return data::spacegroup_type(ds->hall_number);
}

Result<int> SymmetryAnalyzer::spacegroup_number() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->spacegroup_number;
}

Result<int> SymmetryAnalyzer::hall_number() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->hall_number;
}

Result<std::vector<int>> SymmetryAnalyzer::wyckoffs() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->wyckoffs;
}

Result<std::vector<std::string>>
SymmetryAnalyzer::site_symmetry_symbols() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return ds->site_symmetry_symbols;
}

Result<Cell> SymmetryAnalyzer::standardized_cell() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return Cell{ds->std_lattice, ds->std_positions, ds->std_types};
}

Result<Cell> SymmetryAnalyzer::standardized_cell(
    StandardizeOptions options) const {
  return standardize_cell(cell_, options, tol_.symprec, tol_.angle_tolerance);
}

Result<symmetry::Primitive> SymmetryAnalyzer::primitive() const {
  BOOST_LEAF_AUTO(prim, cached_primitive());
  return *prim;
}

Result<Cell> SymmetryAnalyzer::primitive_cell() const {
  BOOST_LEAF_AUTO(prim, cached_primitive());
  return prim->cell;
}

Result<spacegroup::Spacegroup> SymmetryAnalyzer::spacegroup() const {
  BOOST_LEAF_AUTO(sg, cached_spacegroup());
  return *sg;
}

Result<void> SymmetryAnalyzer::warm() const {
  BOOST_LEAF_CHECK(cached_dataset());
  BOOST_LEAF_CHECK(cached_spacegroup()); // also fills cached_primitive()
  BOOST_LEAF_CHECK(cell_operations());
  BOOST_LEAF_CHECK(lattice_symmetry());
  return {};
}

} // namespace spglib::analysis
