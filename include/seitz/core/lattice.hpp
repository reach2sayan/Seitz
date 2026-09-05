#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <cmath>
#include <optional>
#include <utility>

#pragma GCC visibility push(default)

namespace seitz {

// A crystal lattice: the three basis vectors as the columns of a 3x3 matrix.
class Lattice {
public:
  Lattice() = default;
  explicit Lattice(Matrix3d basis) noexcept : basis_{std::move(basis)} {}

  // nullopt for a (near-)singular basis, whose inverse would propagate
  // NaN/inf downstream.
  [[nodiscard]] static std::optional<Lattice>
  from_basis(Matrix3d basis) noexcept {
    if (std::abs(basis.determinant()) < kZeroPrec) {
      return std::nullopt;
    }
    return Lattice{std::move(basis)};
  }

  [[nodiscard]] Matrix3d const &matrix() const noexcept { return basis_; }

  // Cartesian coordinate of a fractional vector, and its inverse.
  [[nodiscard]] Vector3d to_cartesian(Vector3d const &frac) const noexcept {
    return basis_ * frac;
  }
  [[nodiscard]] Vector3d to_fractional(Vector3d const &cart) const noexcept {
    return basis_.inverse() * cart;
  }

  // Metric (Gram) tensor L^T . L, and the geometric volume |det L|.
  [[nodiscard]] Matrix3d metric() const noexcept {
    return basis_.transpose() * basis_;
  }
  [[nodiscard]] double volume() const noexcept {
    return std::abs(basis_.determinant());
  }

  // Change of basis: the lattice whose columns are L . t.
  [[nodiscard]] Lattice transformed(Matrix3d const &t) const noexcept {
    return Lattice{Matrix3d(basis_ * t)};
  }

  // The rotation R with `ideal` = R . *this, from the orthonormal frames of
  // both lattices' first two basis vectors: how the standardized basis is
  // oriented against the Bravais one.
  [[nodiscard]] Matrix3d rigid_rotation_to(Lattice const &ideal) const noexcept;

  // Krivy-Gruber Niggli reduction; `eps` is the tolerance on the metric
  // parameters (typically symprec). Errors e_niggli_failed if it diverges.
  [[nodiscard]] Result<Lattice> niggli(double eps) const;

  // Delaunay reduction, right-handed. Errors e_delaunay_failed on a degenerate
  // cell or a non-unimodular change of basis.
  [[nodiscard]] Result<Lattice> delaunay(double symprec) const;

  // 2D Delaunay reduction in the plane of the two axes other than
  // `unique_axis`, whose column is untouched up to a right-handedness sign
  // flip. The fixed axis is a layer's aperiodic c or a monoclinic unique b.
  [[nodiscard]] Result<Lattice> delaunay_in_plane(int unique_axis,
                                                  double symprec) const;

private:
  Matrix3d basis_{Matrix3d::Identity()};
};

} // namespace seitz

#pragma GCC visibility pop
