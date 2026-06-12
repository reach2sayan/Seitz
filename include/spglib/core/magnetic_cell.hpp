#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/types.hpp>

#include <initializer_list>
#include <utility>
#include <variant>
#include <vector>

namespace spglib {

// Per-site magnetic tensors. The rank is carried by the variant alternative, so
// there is no separate `tensor_rank` field to keep in sync:
//   - collinear     : one scalar magnetic moment per atom.
//   - non-collinear : one 3-vector magnetic moment per atom, stored as an N x 3
//     row-major block (the `Positions` layout) so each moment's three
//     components are contiguous, matching the per-atom symmetry-search access.
using CollinearTensors = std::vector<double>;
using NoncollinearTensors =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using SiteTensors = std::variant<CollinearTensors, NoncollinearTensors>;

// Build noncollinear site tensors from a list of per-atom moments (row i is
// atom i), replacing the brace-init that the std::vector form allowed.
[[nodiscard]] inline NoncollinearTensors
noncollinear_tensors(std::initializer_list<Vector3d> moments) {
  NoncollinearTensors out(static_cast<Index>(moments.size()), 3);
  for (Index i = 0; const Vector3d& m : moments) {
    out.row(i++) = m.transpose();
  }
  return out;
}

// A crystal cell plus per-site magnetic tensors — the input to the magnetic
// symmetry search. Wraps a plain Cell rather than extending it, since site
// tensors are only meaningful for magnetic structures.
class MagneticCell {
public:
  MagneticCell(Cell cell, SiteTensors tensors)
      : cell_(std::move(cell)), tensors_(std::move(tensors)) {}

  [[nodiscard]] Cell const &cell() const noexcept { return cell_; }
  [[nodiscard]] SiteTensors const &tensors() const noexcept { return tensors_; }
  [[nodiscard]] Index size() const noexcept { return cell_.size(); }

  // collinear when the tensors are scalars, non-collinear when 3-vectors.
  [[nodiscard]] SiteTensor rank() const noexcept {
    return std::holds_alternative<CollinearTensors>(tensors_)
               ? SiteTensor::collinear
               : SiteTensor::noncollinear;
  }

  // Scalar / vector moment of atom i. The caller must match the active rank().
  [[nodiscard]] double scalar(Index i) const {
    return std::get<CollinearTensors>(tensors_)[static_cast<std::size_t>(i)];
  }
  [[nodiscard]] Vector3d vector(Index i) const {
    return std::get<NoncollinearTensors>(tensors_).row(i).transpose();
  }

private:
  Cell cell_;
  SiteTensors tensors_;
};

} // namespace spglib
