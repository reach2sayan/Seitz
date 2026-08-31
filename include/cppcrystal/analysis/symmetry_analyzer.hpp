#pragma once

#include <cppcrystal/analysis/detail/lazy.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/dataset.hpp>
#include <cppcrystal/spacegroup/spacegroup.hpp>
#include <cppcrystal/standardize.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>
#include <cppcrystal/symmetry/primitive.hpp>

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
  // Named factory (no overloaded constructors). `hall_number == 0` searches all
  // 230 space groups; a non-zero value fixes the Hall setting.
  [[nodiscard]] static SymmetryAnalyzer
  from_cell(Cell cell, double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt, int hall_number = 0);

  // Analyze `cell` as a layer group with the given aperiodic axis (0/1/2). The
  // axis is stamped onto the owned cell, so dataset()/operations()/etc. all
  // return layer results. (from_cell auto-routes too if the cell already
  // carries an aperiodic axis; this is the explicit convenience form.)
  [[nodiscard]] static SymmetryAnalyzer
  from_layer_cell(Cell cell, int aperiodic_axis,
                  double symprec = kDefaultSymprec,
                  AngleTolerance angle_tolerance = std::nullopt);

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
  [[nodiscard]] Result<SymmetryOperations> operations() const {
    return project<&Dataset::operations>();
  }
  [[nodiscard]] Result<int> spacegroup_number() const {
    return project<&Dataset::spacegroup_number>();
  }
  [[nodiscard]] Result<int> hall_number() const {
    return project<&Dataset::hall_number>();
  }
  [[nodiscard]] Result<std::vector<int>> wyckoffs() const {
    return project<&Dataset::wyckoffs>();
  }
  [[nodiscard]] Result<std::vector<std::string>>
  site_symmetry_symbols() const {
    return project<&Dataset::site_symmetry_symbols>();
  }
  [[nodiscard]] Result<data::SpacegroupType> spacegroup_type() const {
    BOOST_LEAF_AUTO(hall, hall_number());
    return data::spacegroup_type(hall);
  }

  // All space-group operations of the input cell exactly as given, including
  // the centering translations of a non-primitive cell
  // (symmetry::find_symmetry). Distinct from operations(), which are the
  // dataset's operations in the input basis. Cached independently.
  [[nodiscard]] Result<SymmetryOperations> cell_operations() const;

  // The lattice point group: the rotations (in the cell basis) that map the
  // Delaunay-reduced lattice metric onto itself (symmetry::lattice_symmetry).
  [[nodiscard]] Result<PointSymmetry> lattice_symmetry() const;

  // The standardized conventional cell, assembled from the dataset's std_*
  // fields (idealized lattice, fractional positions, atom types).
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // The standardized cell under explicit flags — primitive vs conventional,
  // idealized vs input geometry (standardize_cell). The no-argument overload
  // above is the {} (conventional, idealized) fast path served from the cached
  // dataset; this overload is keyed by `options` and is not memoized.
  [[nodiscard]] Result<Cell> standardized_cell(StandardizeOptions options) const;

  // Pipeline intermediates, cached independently of the full dataset so a
  // caller that only wants the primitive cell or the matched Hall setting does
  // not pay for standardization.
  [[nodiscard]] Result<symmetry::Primitive> primitive() const;
  [[nodiscard]] Result<Cell> primitive_cell() const;
  [[nodiscard]] Result<spacegroup::Spacegroup> spacegroup() const;

  // Force every lazy cache (the dataset plus the independently-cached pipeline
  // intermediates: primitive, matched spacegroup, raw cell operations, lattice
  // symmetry). Returns the first error encountered, or success once all caches
  // are populated. After warm() succeeds this const instance may be shared
  // read-only across threads — all getters then hit a filled cache.
  Result<void> warm() const;

private:
  SymmetryAnalyzer(Cell cell, Tolerance tol, int hall_number)
      : cell_(std::move(cell)), tol_(tol), hall_number_(hall_number) {}

  [[nodiscard]] Result<Dataset const *> cached_dataset() const;
  [[nodiscard]] Result<symmetry::Primitive const *> cached_primitive() const;
  [[nodiscard]] Result<spacegroup::Spacegroup const *>
  cached_spacegroup() const;

  // Copy one field out of the memoized dataset.
  template <auto Member> [[nodiscard]] auto project() const
      -> Result<std::remove_cvref_t<
          decltype(std::declval<Dataset const &>().*Member)>> {
    BOOST_LEAF_AUTO(ds, cached_dataset());
    return ds->*Member;
  }

  Cell cell_;
  Tolerance tol_;
  int hall_number_ = 0;

  detail::Lazy<Dataset> dataset_;
  detail::Lazy<symmetry::Primitive> primitive_;
  detail::Lazy<spacegroup::Spacegroup> spacegroup_;
  detail::Lazy<SymmetryOperations> cell_operations_;
  detail::Lazy<PointSymmetry> lattice_symmetry_;
};

} // namespace cppcrystal::analysis
