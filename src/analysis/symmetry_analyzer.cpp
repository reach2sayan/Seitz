#include <cppcrystal/analysis/symmetry_analyzer.hpp>

#include <utility>

namespace cppcrystal::analysis {

SymmetryAnalyzer SymmetryAnalyzer::from_cell(Cell cell, double symprec,
                                             AngleTolerance angle_tolerance,
                                             int hall_number) {
  return SymmetryAnalyzer{
      std::move(cell),
      Tolerance{.symprec = symprec, .angle_tolerance = angle_tolerance},
      hall_number};
}

SymmetryAnalyzer
SymmetryAnalyzer::from_layer_cell(Cell cell, int aperiodic_axis, double symprec,
                                  AngleTolerance angle_tolerance) {
  cell.set_aperiodic_axis(aperiodic_axis);
  return SymmetryAnalyzer{
      std::move(cell),
      Tolerance{.symprec = symprec, .angle_tolerance = angle_tolerance},
      /*hall_number=*/0};
}

Result<symmetry::Primitive const *> SymmetryAnalyzer::cached_primitive() const {
  return primitive_.get([&] {
    return symmetry::find_primitive(cell_, tol_.symprec, tol_.angle_tolerance);
  });
}

Result<spacegroup::Spacegroup const *>
SymmetryAnalyzer::cached_spacegroup() const {
  return spacegroup_.get([&]() -> Result<spacegroup::Spacegroup> {
    BOOST_LEAF_AUTO(prim, cached_primitive());
    return spacegroup::search_spacegroup(*prim, hall_number_, tol_.symprec,
                                         tol_.angle_tolerance);
  });
}

Result<Dataset const *> SymmetryAnalyzer::cached_dataset() const {
  return dataset_.get([&] {
    return get_dataset(cell_, tol_.symprec, tol_.angle_tolerance, hall_number_);
  });
}

Result<SymmetryOperations> SymmetryAnalyzer::cell_operations() const {
  BOOST_LEAF_AUTO(ops, cell_operations_.get([&] {
    return symmetry::find_symmetry(cell_, tol_.symprec, tol_.angle_tolerance);
  }));
  return *ops;
}

Result<PointSymmetry> SymmetryAnalyzer::lattice_symmetry() const {
  BOOST_LEAF_AUTO(ps, lattice_symmetry_.get([&] {
    return symmetry::lattice_symmetry(cell_, tol_.symprec,
                                      tol_.angle_tolerance);
  }));
  return *ps;
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

} // namespace cppcrystal::analysis
