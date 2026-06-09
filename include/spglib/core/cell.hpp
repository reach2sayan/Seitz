#pragma once

#include <spglib/core/types.hpp>
#include <spglib/math/integer_matrix.hpp>

#include <cmath>
#include <optional>
#include <utility>

namespace spglib {

// Rank of per-site tensors for magnetic structures (spglib SiteTensorType).
enum class SiteTensor { none = -1, collinear = 0, noncollinear = 1 };

// A crystal cell: lattice with columns equal to the Cartesian basis vectors,
// fractional atomic positions (row i is atom i), and integer atom types.
class Cell {
public:
  Cell() = default;
  Cell(Matrix3d lattice, Positions positions, Types types,
       std::optional<int> aperiodic_axis = std::nullopt)
      : lattice_(std::move(lattice)), positions_(std::move(positions)),
        types_(std::move(types)), aperiodic_axis_(aperiodic_axis) {}

  [[nodiscard]] Index size() const noexcept {
    return static_cast<Index>(types_.size());
  }

  [[nodiscard]] Matrix3d const &lattice() const noexcept { return lattice_; }
  [[nodiscard]] Positions const &positions() const noexcept {
    return positions_;
  }
  [[nodiscard]] Types const &types() const noexcept { return types_; }

  // std::nullopt for a 3D space-group search; the aperiodic axis index (0..2)
  // for a layer group.
  [[nodiscard]] std::optional<int> aperiodic_axis() const noexcept {
    return aperiodic_axis_;
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
    aperiodic_axis_ = axis;
  }

private:
  Matrix3d lattice_{Matrix3d::Identity()};
  Positions positions_{};
  Types types_{};
  std::optional<int> aperiodic_axis_;
};

} // namespace spglib
