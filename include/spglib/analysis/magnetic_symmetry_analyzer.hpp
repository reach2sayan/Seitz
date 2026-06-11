#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/magnetic_cell.hpp>
#include <spglib/core/magnetic_symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/data/msg_database.hpp>
#include <spglib/magnetic_dataset.hpp>

#include <optional>
#include <vector>

namespace spglib::analysis {

// The magnetic counterpart of SymmetryAnalyzer: a persistent view over a
// MagneticCell + tolerances that lazily computes and memoizes the magnetic
// space-group analysis by delegating to get_magnetic_dataset. Same design as
// SymmetryAnalyzer — immutable after construction, value-owned inputs, the
// success value cached (boost::leaf::result is move-only), not thread-safe.
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
  [[nodiscard]] Result<int> uni_number() const;
  [[nodiscard]] Result<int> hall_number() const;
  [[nodiscard]] Result<data::MagneticSpacegroupType> spacegroup_type() const;
  [[nodiscard]] Result<std::vector<int>> equivalent_atoms() const;

  // Standardized magnetic cell assembled from the dataset's std_* fields
  // (idealized lattice, positions, types, and the rotated site tensors).
  [[nodiscard]] Result<MagneticCell> standardized_cell() const;

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
};

} // namespace spglib::analysis
