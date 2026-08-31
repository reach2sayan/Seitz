#pragma once

// Eigen/Dense (not just Eigen/Core) so that determinant(), inverse(), and the
// other decompositions are defined for our matrix aliases in every translation
// unit that uses them. Including only Eigen/Core leaves those as undefined
// references that may or may not be satisfied by another TU's instantiation.
#include <Eigen/Dense>

#include <span>
#include <vector>

namespace cppcrystal {

using Index = Eigen::Index;

using Vector3d = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using Matrix3d = Eigen::Matrix3d; // columns are the (Cartesian) basis vectors
using Matrix3i = Eigen::Matrix3i; // integer rotation expressed in the lat basis

// N x 3 block of fractional (or Cartesian) coordinates; row i is point/atom i.
// Row-major so each point's three components are contiguous, matching the
// per-atom access patterns that dominate the symmetry search.
using Positions = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using Types = std::vector<int>;

// Index-typed access to a Types vector. Eigen::Index is signed and
// std::vector indexes by size_t — this is the single place that cast lives.
[[nodiscard]] inline int &type_at(Types &types, Index i) noexcept {
  return types[static_cast<std::size_t>(i)];
}
[[nodiscard]] inline int type_at(Types const &types, Index i) noexcept {
  return types[static_cast<std::size_t>(i)];
}

// Pack contiguous Vector3d rows into an N x 3 Positions block in one copy;
// Positions is RowMajor and Vector3d rows are laid out exactly the same way.
[[nodiscard]] inline Positions to_positions(std::span<Vector3d const> rows) {
  static_assert(sizeof(Vector3d) == 3 * sizeof(double));
  if (rows.empty()) {
    return Positions{0, 3};
  }
  return Eigen::Map<Positions const>(rows.data()->data(),
                                     static_cast<Index>(rows.size()), 3);
}

} // namespace cppcrystal
