#pragma once

// Hand-written Eigen views over the generated hall_symbol.c data tables
// (hall_generators.hpp). The raw tables stay plain constexpr std::array so they
// live in .rodata at compile time; these accessors are the single place that
// bridges them into Eigen. They pin the row-major layout of the flattened
// storage so call sites cannot reintroduce a silent transpose (a bare
// Eigen::Map<Matrix...> over this buffer would read it column-major).

#include <spglib/core/types.hpp>
#include <spglib/data/hall_generators.hpp>

#include <cstddef>

namespace spglib::data {

// A VSpU set is a 3x9 pseudo-inverse (V.Sigma+.U^T) stored as three contiguous
// rows of 9 -> row-major 3x9. shift = VSpU . dw.
using VSpUMatrix = Eigen::Matrix<double, 3, 9, Eigen::RowMajor>;

// The dw column (per-generator translation differences) that VSpU multiplies.
using DwVector = Eigen::Matrix<double, 9, 1>;

[[nodiscard]] inline Eigen::Map<VSpUMatrix const>
vspu_matrix(VSpUSet const &v) {
  // std::array<std::array<double,9>,3> is contiguous, so v[0].data() is the
  // start of all 27 row-major values.
  return Eigen::Map<VSpUMatrix const>(v[0].data());
}

// One generator: a flattened 3x3 rotation, row-major.
using GeneratorMatrix = Eigen::Matrix<int, 3, 3, Eigen::RowMajor>;

[[nodiscard]] inline Eigen::Map<GeneratorMatrix const>
generator_matrix(GeneratorSet const &g, std::size_t i) {
  return Eigen::Map<GeneratorMatrix const>(g[i].data());
}

} // namespace spglib::data
