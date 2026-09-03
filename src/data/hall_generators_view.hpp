#pragma once

#include "data/hall_generators.hpp"
#include <cppcrystal/core/types.hpp>

#include <cstddef>

namespace cppcrystal::data {

// A VSpU set is a 3x9 pseudo-inverse (V.Sigma+.U^T) stored as three contiguous
// rows of 9 -> row-major 3x9. shift = VSpU . dw.
using VSpUMatrix = Eigen::Matrix<double, 3, 9, Eigen::RowMajor>;

// The dw column (per-generator translation differences) that VSpU multiplies.
using DwVector = Eigen::Matrix<double, 9, 1>;

[[nodiscard]] inline Eigen::Map<VSpUMatrix const>
vspu_matrix(VSpUSet const &v) {
  return Eigen::Map<VSpUMatrix const>(v[0].data());
}

// One generator: a flattened 3x3 rotation, row-major.
using GeneratorMatrix = Eigen::Matrix<int, 3, 3, Eigen::RowMajor>;

[[nodiscard]] inline Eigen::Map<GeneratorMatrix const>
generator_matrix(GeneratorSet const &g, std::size_t i) {
  return Eigen::Map<GeneratorMatrix const>(g[i].data());
}

} // namespace cppcrystal::data
