#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/data/msg_database.hpp>
#include <spglib/magnetic_dataset.hpp>
#include <spglib/spin/spin.hpp>

#include <optional>
#include <vector>

namespace spglib::analysis {

// The magnetic counterpart of SymmetryAnalyzer: a persistent view over a
// MagneticCell + tolerances that lazily computes and memoizes the magnetic
// space-group analysis via get_magnetic_dataset. Same design as SymmetryAnalyzer
// — immutable after construction, value-owned inputs, the success value cached.
// The per-instance caches are not thread-safe; call warm() once on a single
// thread to make the const instance shareable read-only across threads.
class MagneticSymmetryAnalyzer {
public:
  // `is_axial` selects axial-vector transformation of the rank-1 tensors;
  // `mag_symprec` (std::nullopt -> symprec) is the moment tolerance.
  [[nodiscard]] static MagneticSymmetryAnalyzer
  from_cell(MagneticCell cell, bool is_axial,
            double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt,
            std::optional<double> mag_symprec = std::nullopt);

  [[nodiscard]] MagneticCell const &cell() const noexcept { return cell_; }
  [[nodiscard]] bool is_axial() const noexcept { return is_axial_; }
  [[nodiscard]] double symprec() const noexcept { return tol_.symprec; }
  [[nodiscard]] AngleTolerance angle_tolerance() const noexcept {
    return tol_.angle_tolerance;
  }

  [[nodiscard]] Result<MagneticDataset> dataset() const;

  [[nodiscard]] Result<MagneticSymmetryOperations> operations() const;

  // The full magnetic symmetry search of the input cell: the magnetic
  // operations, each atom's magnetic-orbit representative, the atom permutations
  // under every operation, and the primitive lattice implied by the magnetic
  // pure translations. Derives the spatial symmetry (symmetry::find_symmetry of
  // the underlying cell) and runs spin::operations_with_site_tensors with time
  // reversal. Distinct from operations(), which is just the dataset's magnetic
  // operations. Cached independently.
  [[nodiscard]] Result<spin::MagneticSymmetrySearch> symmetry_search() const;
  [[nodiscard]] Result<int> uni_number() const;
  [[nodiscard]] Result<int> hall_number() const;
  [[nodiscard]] Result<data::MagneticSpacegroupType> spacegroup_type() const;
  [[nodiscard]] Result<std::vector<int>> equivalent_atoms() const;

  // Standardized magnetic cell assembled from the dataset's std_* fields
  // (idealized lattice, positions, types, and the rotated site tensors).
  [[nodiscard]] Result<MagneticCell> standardized_cell() const;

  // Force both lazy caches (the magnetic dataset and the independently-cached
  // symmetry search). After warm() succeeds this const instance may be shared
  // read-only across threads.
  Result<void> warm() const;

private:
  MagneticSymmetryAnalyzer(MagneticCell cell, bool is_axial, Tolerance tol,
                           std::optional<double> mag_symprec)
      : cell_(std::move(cell)), is_axial_(is_axial), tol_(tol),
        mag_symprec_(mag_symprec) {}

  [[nodiscard]] Result<MagneticDataset const *> cached_dataset() const;

  MagneticCell cell_;
  bool is_axial_ = false;
  Tolerance tol_;
  std::optional<double> mag_symprec_;

  mutable std::optional<MagneticDataset> dataset_;
  mutable std::optional<spin::MagneticSymmetrySearch> symmetry_search_;
};

} // namespace spglib::analysis
