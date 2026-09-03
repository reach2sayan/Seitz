#pragma once

#include <cppcrystal/analysis/detail/lazy.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/dataset.hpp>
#include <cppcrystal/standardize.hpp>

#include <optional>
#include <string>
#include <type_traits>
#include <vector>

namespace cppcrystal::analysis {

// A persistent, stateful view over a Cell + tolerances that lazily computes and
// memoizes the symmetry analysis: owns the inputs once and caches each pipeline
// stage (detail::Lazy) so repeated queries do not recompute. Immutable after
// construction (no setters); to analyze at a different tolerance, build a new
// analyzer.
//
// Thread-safety: the per-instance caches are NOT race-free — concurrent
// first-calls race. To share one instance read-only across threads, call
// warm() once on a single thread first; afterwards every getter is served from
// a populated cache and does no writing. (The shared global tables are primed
// separately by cppcrystal::warmup().)
class SymmetryAnalyzer {
public:
  // Named factory (no overloaded constructors). An unset `setting` searches
  // every Hall setting of the cell's family; a set one fixes it.
  [[nodiscard]] static SymmetryAnalyzer
  from_cell(Cell cell, Tolerance tol = {},
            std::optional<HallNumber> setting = std::nullopt);

  [[nodiscard]] Cell const &cell() const noexcept { return cell_; }
  [[nodiscard]] double symprec() const noexcept { return tol_.symprec; }
  [[nodiscard]] AngleTolerance angle_tolerance() const noexcept {
    return tol_.angle_tolerance;
  }

  // The full space-group dataset of the input cell (memoized). Every getter
  // below is a projection of this and shares the one computation.
  [[nodiscard]] Result<Dataset> dataset() const {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return *ds;
  }

  // Projections of the dataset.
  [[nodiscard]] Result<Operations> operations() const {
    return project<&Dataset::operations>();
  }
  [[nodiscard]] Result<int> spacegroup_number() const {
    return project<&Dataset::spacegroup_number>();
  }
  [[nodiscard]] Result<HallNumber> hall() const {
    return project<&Dataset::hall>();
  }
  [[nodiscard]] Result<std::vector<int>> wyckoffs() const {
    return project<&Dataset::wyckoffs>();
  }
  [[nodiscard]] Result<std::vector<std::string>>
  site_symmetry_symbols() const {
    return project<&Dataset::site_symmetry_symbols>();
  }
  [[nodiscard]] Result<data::SpacegroupType> spacegroup_type() const {
    BOOST_LEAF_AUTO(setting, hall());
    return data::spacegroup_type(setting);
  }

  // All space-group operations of the input cell exactly as given, including
  // the centering translations of a non-primitive cell
  // (symmetry::find_symmetry). Distinct from operations(), which are the
  // dataset's operations in the input basis. Cached independently.
  [[nodiscard]] Result<Operations> cell_operations() const;

  // The lattice point group: the rotations (in the cell basis) that map the
  // Delaunay-reduced lattice metric onto itself (symmetry::lattice_symmetry).
  [[nodiscard]] Result<PointSymmetry> lattice_symmetry() const;

  // The standardized conventional cell, assembled from the dataset's std_*
  // fields (idealized lattice, fractional positions, atom types).
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // The standardized cell in an explicit setting — primitive vs conventional,
  // idealized vs input geometry (standardize_cell). The no-argument overload
  // above is the (conventional, idealized) fast path served from the cached
  // dataset; this one is keyed by its arguments and is not memoized.
  [[nodiscard]] Result<Cell> standardized_cell(CellSetting setting,
                                               Idealize idealize) const;

  // The primitive cell, cached independently of the full dataset so a caller
  // that only wants it does not pay for standardization. The Primitive and
  // Spacegroup intermediates behind it are private to the pipeline.
  [[nodiscard]] Result<Cell> primitive_cell() const;

  // Force every lazy cache (the dataset plus the independently-cached
  // primitive cell, raw cell operations and lattice symmetry). Returns the
  // first error encountered, or success once all caches
  // are populated. After warm() succeeds this const instance may be shared
  // read-only across threads — all getters then hit a filled cache.
  Result<void> warm() const;

private:
  SymmetryAnalyzer(Cell cell, Tolerance tol, std::optional<HallNumber> setting)
      : cell_(std::move(cell)), tol_(tol), setting_(setting) {}

  [[nodiscard]] Result<Dataset const *> cached_dataset() const;

  // Copy one field out of the memoized dataset.
  template <auto Member> [[nodiscard]] auto project() const
      -> Result<std::remove_cvref_t<
          decltype(std::declval<Dataset const &>().*Member)>> {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return ds->*Member;
  }

  Cell cell_;
  Tolerance tol_;
  std::optional<HallNumber> setting_;

  detail::Lazy<Dataset> dataset_;
  detail::Lazy<Cell> primitive_cell_;
  detail::Lazy<Operations> cell_operations_;
  detail::Lazy<PointSymmetry> lattice_symmetry_;
};

} // namespace cppcrystal::analysis
