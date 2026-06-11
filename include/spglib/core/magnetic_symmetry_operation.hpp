#pragma once

#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/types.hpp>

#include <vector>

namespace spglib {

// A magnetic space-group operation: an ordinary space-group operation paired
// with a time-reversal flag. `time_reversal == true` is an anti-operation (the
// "primed" element that also reverses time / flips magnetic moments); `false`
// is an ordinary unitary operation. Replaces spglib's parallel
// MagneticSymmetry.rot[] / .trans[] / .timerev[] arrays with one value type,
// collected in a std::vector (see MagneticSymmetryOperations). The spglib C
// convention encodes the flag as the int timerev (1 = anti, 0 = ordinary); we
// keep it a bool — there is no third state.
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

} // namespace spglib
