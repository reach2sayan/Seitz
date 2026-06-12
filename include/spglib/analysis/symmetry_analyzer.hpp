#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/dataset.hpp>
#include <spglib/spacegroup/spacegroup.hpp>
#include <spglib/standardize.hpp>
#include <spglib/symmetry/find_symmetry.hpp>
#include <spglib/symmetry/primitive.hpp>

#include <optional>
#include <string>
#include <vector>

namespace spglib::analysis {

// A persistent, stateful view over a Cell + tolerances that lazily computes and
// memoizes the symmetry analysis — the modern-C++ analogue of pymatgen's
// SpacegroupAnalyzer. The functional core (get_dataset, find_primitive,
// search_spacegroup) is left untouched; this facade owns the inputs once and
// caches each pipeline stage so that repeated queries do not recompute.
//
// The analyzer is immutable after construction: inputs are owned by value and
// there are no setters. To analyze at a different tolerance, build a new
// analyzer. Memoization makes the const getters logically const; the caches are
// `mutable`.
//
// Thread-safety: the shared global tables this delegates to (the Hall operation
// database etc.) are primed once by spglib::warmup() and are then race-free for
// concurrent reads. The per-instance caches are NOT — concurrent first-calls race
// on the `mutable` optionals. To share one instance read-only across threads,
// call warm() once on a single thread first; afterwards every getter is served
// from a populated cache and does no writing.
//
// boost::leaf::result is move-only, so the caches store the success *value*
// (std::optional<T>) rather than the Result. Each getter returns a freshly
// constructed Result; on error nothing is cached and the next call re-runs
// (errors are exceptional, so this costs nothing on the happy path).
class SymmetryAnalyzer {
public:
  // PyXtal-style named factory (no overloaded constructors). `hall_number == 0`
  // searches all 230 space groups; a non-zero value fixes the Hall setting.
  [[nodiscard]] static SymmetryAnalyzer
  from_cell(Cell cell, double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt, int hall_number = 0);

  // Analyze `cell` as a layer group with the given aperiodic axis (0/1/2). The
  // axis is stamped onto the owned cell, so dataset()/operations()/etc. all
  // return layer results. (from_cell auto-routes too if the cell already carries
  // an aperiodic axis; this is the explicit convenience form.)
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
  [[nodiscard]] Result<Dataset> dataset() const;

  // Projections of the dataset.
  [[nodiscard]] Result<SymmetryOperations> operations() const;

  // All space-group operations of the input cell exactly as given, including the
  // centering translations of a non-primitive cell (symmetry::find_symmetry —
  // the raw sym_get_operation result). Distinct from operations(), which are the
  // dataset's operations in the input basis. Cached independently.
  [[nodiscard]] Result<SymmetryOperations> cell_operations() const;

  // The lattice point group: the rotations (in the cell basis) that map the
  // Delaunay-reduced lattice metric onto itself (symmetry::lattice_symmetry).
  [[nodiscard]] Result<PointSymmetry> lattice_symmetry() const;
  [[nodiscard]] Result<data::SpacegroupType> spacegroup_type() const;
  [[nodiscard]] Result<int> spacegroup_number() const;
  [[nodiscard]] Result<int> hall_number() const;
  [[nodiscard]] Result<std::vector<int>> wyckoffs() const;
  [[nodiscard]] Result<std::vector<std::string>> site_symmetry_symbols() const;

  // The standardized conventional cell, assembled from the dataset's std_*
  // fields (idealized lattice, fractional positions, atom types).
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // The standardized cell under explicit flags — primitive vs conventional,
  // idealized vs input geometry (standardize_cell). The no-argument overload
  // above is the {} (conventional, idealized) fast path served from the cached
  // dataset; this overload is keyed by `options` and is not memoized.
  [[nodiscard]] Result<Cell> standardized_cell(StandardizeOptions options) const;

  // Pipeline intermediates, cached independently of the full dataset so a caller
  // that only wants the primitive cell or the matched Hall setting does not pay
  // for standardization.
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

  // Ensure-and-return helpers: populate the cache on first use and hand back a
  // pointer into it (or propagate the error). Returning a pointer lets the
  // projections copy out only the field they need rather than the whole value.
  [[nodiscard]] Result<Dataset const *> cached_dataset() const;
  [[nodiscard]] Result<symmetry::Primitive const *> cached_primitive() const;
  [[nodiscard]] Result<spacegroup::Spacegroup const *> cached_spacegroup() const;

  Cell cell_;
  Tolerance tol_;
  int hall_number_ = 0;

  mutable std::optional<Dataset> dataset_;
  mutable std::optional<symmetry::Primitive> primitive_;
  mutable std::optional<spacegroup::Spacegroup> spacegroup_;
  mutable std::optional<SymmetryOperations> cell_operations_;
  mutable std::optional<PointSymmetry> lattice_symmetry_;
};

} // namespace spglib::analysis
