#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/dataset.hpp>
#include <spglib/spacegroup/spacegroup.hpp>
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
// `mutable`. NOT thread-safe — concurrent first-calls race on the cache; treat
// an instance as single-threaded until warmed.
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
  [[nodiscard]] Result<data::SpacegroupType> spacegroup_type() const;
  [[nodiscard]] Result<int> spacegroup_number() const;
  [[nodiscard]] Result<int> hall_number() const;
  [[nodiscard]] Result<std::vector<int>> wyckoffs() const;
  [[nodiscard]] Result<std::vector<std::string>> site_symmetry_symbols() const;

  // The standardized conventional cell, assembled from the dataset's std_*
  // fields (idealized lattice, fractional positions, atom types).
  [[nodiscard]] Result<Cell> standardized_cell() const;

  // Pipeline intermediates, cached independently of the full dataset so a caller
  // that only wants the primitive cell or the matched Hall setting does not pay
  // for standardization.
  [[nodiscard]] Result<symmetry::Primitive> primitive() const;
  [[nodiscard]] Result<Cell> primitive_cell() const;
  [[nodiscard]] Result<spacegroup::Spacegroup> spacegroup() const;

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
};

} // namespace spglib::analysis
