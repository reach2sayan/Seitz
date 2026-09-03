#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>

#include "core/family_dispatch.hpp"
#include "symmetry/primitive.hpp"
#include "symmetry/search.hpp"

#include <utility>

namespace cppcrystal::analysis {

SymmetryAnalyzer SymmetryAnalyzer::from_cell(Cell cell, Tolerance tol,
                                             int hall_number) {
  return SymmetryAnalyzer{std::move(cell), tol, hall_number};
}

Result<Dataset const *> SymmetryAnalyzer::cached_dataset() const {
  return dataset_.get([&] { return get_dataset(cell_, tol_, hall_number_); });
}

Result<Operations> SymmetryAnalyzer::cell_operations() const {
  BOOST_LEAF_AUTO(ops, cell_operations_.get([&] {
    return dispatch_family(cell_.periodicity(), [&]<GroupFamily F>() {
      return symmetry::SymmetrySearch<F>{cell_, tol_}.operations();
    });
  }));
  return *ops;
}

Result<PointSymmetry> SymmetryAnalyzer::lattice_symmetry() const {
  BOOST_LEAF_AUTO(ps, lattice_symmetry_.get([&] {
    return dispatch_family(cell_.periodicity(), [&]<GroupFamily F>() {
      return symmetry::SymmetrySearch<F>{cell_, tol_}.lattice_symmetry();
    });
  }));
  return *ps;
}

Result<Cell> SymmetryAnalyzer::standardized_cell() const {
  BOOST_LEAF_AUTO(ds, cached_dataset());
  return Cell{Lattice{ds->std_lattice}, ds->std_positions, ds->std_types};
}

Result<Cell> SymmetryAnalyzer::standardized_cell(CellSetting setting,
                                                 Idealize idealize) const {
  return standardize_cell(cell_, setting, idealize, tol_);
}

Result<Cell> SymmetryAnalyzer::primitive_cell() const {
  BOOST_LEAF_AUTO(cell, primitive_cell_.get([&]() -> Result<Cell> {
    return dispatch_family(
        cell_.periodicity(), [&]<GroupFamily F>() -> Result<Cell> {
          symmetry::PrimitiveFinder<F> const finder(cell_, tol_);
          BOOST_LEAF_AUTO(prim, finder.find());
          return prim.cell;
        });
  }));
  return *cell;
}

Result<void> SymmetryAnalyzer::warm() const {
  BOOST_LEAF_CHECK(cached_dataset());
  BOOST_LEAF_CHECK(primitive_cell());
  BOOST_LEAF_CHECK(cell_operations());
  BOOST_LEAF_CHECK(lattice_symmetry());
  return {};
}

} // namespace cppcrystal::analysis
