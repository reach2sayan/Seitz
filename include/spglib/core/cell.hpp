#pragma once

#include <spglib/core/periodicity.hpp>
#include <spglib/core/types.hpp>
#include <spglib/math/integer_matrix.hpp>

#include <cmath>
#include <optional>
#include <utility>

namespace spglib {

// Rank of per-site tensors for magnetic structures.
enum class SiteTensor { none = -1, collinear = 0, noncollinear = 1 };

// A crystal cell: lattice with columns equal to the Cartesian basis vectors,
// fractional atomic positions (row i is atom i), and integer atom types.
class Cell {
public:
  Cell() = default;
  // The legacy constructor: a 3D cell (default) or a single aperiodic axis (a
  // layer group's c). Retained so existing call sites are unchanged; it sets the
  // per-axis periodicity accordingly.
  Cell(Matrix3d lattice, Positions positions, Types types,
       std::optional<int> aperiodic_axis = std::nullopt)
      : lattice_(std::move(lattice)), positions_(std::move(positions)),
        types_(std::move(types)),
        periodicity_(periodicity_from_aperiodic_axis(aperiodic_axis)) {}
  // The general constructor: an explicit per-axis periodicity, the only form
  // that can describe a rod (two aperiodic axes) or a 0D cluster.
  Cell(Matrix3d lattice, Positions positions, Types types,
       CellPeriodicity periodicity)
      : lattice_(std::move(lattice)), positions_(std::move(positions)),
        types_(std::move(types)), periodicity_(periodicity) {}

  [[nodiscard]] Index size() const noexcept {
    return static_cast<Index>(types_.size());
  }

  [[nodiscard]] Matrix3d const &lattice() const noexcept { return lattice_; }
  [[nodiscard]] Positions const &positions() const noexcept {
    return positions_;
  }
  [[nodiscard]] Types const &types() const noexcept { return types_; }

  // The full per-axis periodicity (the canonical descriptor).
  [[nodiscard]] CellPeriodicity const &periodicity() const noexcept {
    return periodicity_;
  }

  // The single aperiodic axis if there is exactly one (a layer group's c),
  // std::nullopt otherwise — i.e. for a 3D cell (all periodic) AND for the
  // rod/cluster cases with more than one aperiodic axis, which this legacy view
  // cannot represent (use periodicity() there). Backward-compatible with the
  // previous std::optional<int> field for space-group / layer code.
  [[nodiscard]] std::optional<int> aperiodic_axis() const noexcept {
    return single_aperiodic_axis(periodicity_);
  }

  // Fractional position / type of atom i.
  [[nodiscard]] Vector3d position(Index i) const noexcept {
    return positions_.row(i).transpose();
  }
  [[nodiscard]] int type(Index i) const noexcept {
    return types_[static_cast<std::size_t>(i)];
  }

  // Cartesian coordinate of a fractional vector: L . frac.
  [[nodiscard]] Vector3d to_cartesian(Vector3d const &frac) const noexcept {
    return lattice_ * frac;
  }

  // Geometric cell volume |det(L)|.
  [[nodiscard]] double volume() const noexcept {
    return std::abs(lattice_.determinant());
  }

  // Metric (Gram) tensor L^T . L.
  [[nodiscard]] Matrix3d metric() const noexcept {
    return math::metric_tensor(lattice_);
  }

  // Mutable access for cell builders (primitive/standardized construction).
  [[nodiscard]] Matrix3d &lattice() noexcept { return lattice_; }
  [[nodiscard]] Positions &positions() noexcept { return positions_; }
  [[nodiscard]] Types &types() noexcept { return types_; }
  void set_aperiodic_axis(std::optional<int> axis) noexcept {
    periodicity_ = periodicity_from_aperiodic_axis(axis);
  }
  void set_periodicity(CellPeriodicity periodicity) noexcept {
    periodicity_ = periodicity;
  }

private:
  Matrix3d lattice_{Matrix3d::Identity()};
  Positions positions_{};
  Types types_{};
  CellPeriodicity periodicity_{all_periodic()};
};

} // namespace spglib
