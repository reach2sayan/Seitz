#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/types.hpp>

#include <initializer_list>
#include <utility>
#include <variant>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz {

// How a rank-1 site tensor transforms under an improper operation:
// (a) axial : picks up the |det R| factor,
// (b) polar : one does not.
enum class TensorKind { polar, axial };

// Rank of per-site tensors for magnetic structures.
enum class SiteTensor { none = -1, collinear = 0, noncollinear = 1 };

// Per-site magnetic tensors; the rank is the variant alternative, so no
// `tensor_rank` field can fall out of sync:
//   collinear     : one scalar moment per atom.
//   non-collinear : one 3-vector moment per atom, as an N x 3 row-major block
//                   (the `Positions` layout), contiguous per atom.
using CollinearTensors = std::vector<double>;
using NoncollinearTensors =
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using SiteTensors = std::variant<CollinearTensors, NoncollinearTensors>;

// Noncollinear site tensors from per-atom moments, row i = atom i.
[[nodiscard]] inline NoncollinearTensors
noncollinear_tensors(std::initializer_list<Vector3d> moments) {
  NoncollinearTensors out(static_cast<Index>(moments.size()), 3);
  for (Index i = 0; const Vector3d &m : moments) {
    out.row(i++) = m.transpose();
  }
  return out;
}

// A Cell plus per-site magnetic tensors: the input to the magnetic symmetry
// search. Wraps a Cell rather than extending it -- site tensors mean nothing
// for a non-magnetic structure.
class MagneticCell {
public:
  MagneticCell(Cell cell, SiteTensors tensors,
               TensorKind kind = TensorKind::polar)
      : cell_(std::move(cell)), tensors_(std::move(tensors)), kind_(kind) {}

  [[nodiscard]] Cell const &cell() const noexcept { return cell_; }
  [[nodiscard]] SiteTensors const &tensors() const noexcept { return tensors_; }
  [[nodiscard]] TensorKind kind() const noexcept { return kind_; }
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
  TensorKind kind_;
};

} // namespace seitz

#pragma GCC visibility pop
