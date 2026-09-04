#include <cppcrystal/core/lattice.hpp>

namespace cppcrystal {

namespace {

// Orthonormal frame from the first two basis vectors
// (e1 = a^, e3 = (a x b)^, e2 = e3 x e1); columns are the frame vectors.
// e2 needs no normalising: e3 and e1 are unit and orthogonal, so their cross
// product is already unit.
[[nodiscard]] Matrix3d orthonormal_basis(Matrix3d const &lattice) noexcept {
  Vector3d const a = lattice.col(0);
  Vector3d const b = lattice.col(1);
  Vector3d const e1 = a.normalized();
  Vector3d const e3 = a.cross(b).normalized();
  Vector3d const e2 = e3.cross(e1);

  Matrix3d basis;
  basis << e1, e2, e3;
  return basis;
}

} // namespace

Matrix3d Lattice::rigid_rotation_to(Lattice const &ideal) const noexcept {
  // The frames are orthonormal, so the inverse is the transpose -- exact, and
  // without running a 3x3 cofactor inversion to arrive at the same thing.
  return orthonormal_basis(ideal.basis_) * orthonormal_basis(basis_).transpose();
}

} // namespace cppcrystal
