#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <vector>

namespace cppcrystal {

// A magnetic space-group operation: an ordinary space-group operation paired
// with a time-reversal flag. `time_reversal == true` is an anti-operation (the
// "primed" element that also reverses time / flips magnetic moments); `false`
// is an ordinary unitary operation.
struct MagneticSymmetryOperation {
  Matrix3i rotation{Matrix3i::Identity()};
  Vector3d translation{Vector3d::Zero()};
  bool time_reversal = false;

  // The underlying space-group operation, dropping the time-reversal flag.
  [[nodiscard]] SymmetryOperation spatial() const noexcept {
    return {rotation, translation};
  }
};

using MagneticSymmetryOperations = std::vector<MagneticSymmetryOperation>;

} // namespace cppcrystal
