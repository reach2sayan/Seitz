#pragma once

#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#pragma GCC visibility push(default)

namespace seitz {

// A magnetic space-group operation: an ordinary space-group operation paired
// with a time-reversal flag. `time_reversal == true` is an anti-operation (the
// "primed" element that also reverses time / flips magnetic moments); `false`
// is an ordinary unitary operation.
struct MagneticSymmetryOperation {
  SymmetryOperation spatial{};
  bool time_reversal = false;
};

} // namespace seitz

#pragma GCC visibility pop
