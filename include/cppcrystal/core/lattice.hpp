#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include <cmath>
#include <optional>
#include <utility>

namespace cppcrystal {

// A crystal lattice: the three basis vectors as the columns of a 3x3 matrix.
//
// matrix() is the ONE Eigen escape hatch — every other operation the pipeline
// performs on "a lattice" (metric, volume, coordinate conversion, reduction,
// change of basis) is a member here, so a bare Matrix3d never has to carry the
// implicit "these columns are basis vectors" contract on its own.
class Lattice {
public:
  Lattice() = default;
  explicit Lattice(Matrix3d basis) noexcept : basis_(std::move(basis)) {}

  // Checked construction: std::nullopt for a (near-)singular basis, whose
  // inverse would otherwise propagate NaN/inf through the whole pipeline.
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

  // Rigid rotation R with `ideal` = R . *this, built from the orthonormal
  // frames of the two lattices' first two basis vectors. Used to record how the
  // standardized basis is oriented relative to the bravais one.
  [[nodiscard]] Matrix3d rigid_rotation_to(Lattice const &ideal) const noexcept;

  // Krivy-Gruber Niggli reduction. `eps` is the tolerance on the metric
  // parameters (typically symprec). Errors with e_niggli_failed if it does not
  // converge.
  [[nodiscard]] Result<Lattice> niggli(double eps) const;

  // Delaunay reduction; the result is right-handed. Errors with
  // e_delaunay_failed when the cell is degenerate or the change of basis is not
  // unimodular.
  [[nodiscard]] Result<Lattice> delaunay(double symprec) const;

  // 2D Delaunay reduction within the plane spanned by the two axes other than
  // `unique_axis`, whose column is left untouched (up to a sign flip that keeps
  // the result right-handed). The extra argument distinguishes this from the 3D
  // overload above.
  [[nodiscard]] Result<Lattice> delaunay(int unique_axis, double symprec) const;

private:
  Matrix3d basis_{Matrix3d::Identity()};
};

} // namespace cppcrystal
