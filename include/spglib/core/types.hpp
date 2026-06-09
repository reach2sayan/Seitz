#pragma once

// Eigen/Dense (not just Eigen/Core) so that determinant(), inverse(), and the
// other decompositions are defined for our matrix aliases in every translation
// unit that uses them. Including only Eigen/Core leaves those as undefined
// references that may or may not be satisfied by another TU's instantiation.
#include <Eigen/Dense>

#include <vector>

namespace spglib {

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

} // namespace spglib
