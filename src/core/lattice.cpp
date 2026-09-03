#include <cppcrystal/core/lattice.hpp>

namespace cppcrystal {

namespace {

// Orthonormal frame from the first two basis vectors
// (e1 = a^, e3 = (a x b)^, e2 = (e3 x a)^); columns are the frame vectors.
[[nodiscard]] Matrix3d orthonormal_basis(Matrix3d const &lattice) noexcept {
  Vector3d const a = lattice.col(0);
  Vector3d const b = lattice.col(1);
  Vector3d const e1 = a.normalized();
  Vector3d const e3 = a.cross(b).normalized();
  Vector3d const e2 = e3.cross(a).normalized();

  Matrix3d basis;
  basis << e1, e2, e3;
  return basis;
}

} // namespace

Matrix3d Lattice::rigid_rotation_to(Lattice const &ideal) const noexcept {
  return orthonormal_basis(ideal.basis_) * orthonormal_basis(basis_).inverse();
}

} // namespace cppcrystal
