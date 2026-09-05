#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#pragma GCC visibility push(default)

namespace cppcrystal {

// A magnetic space-group operation: an ordinary space-group operation paired
// with a time-reversal flag. `time_reversal == true` is an anti-operation (the
// "primed" element that also reverses time / flips magnetic moments); `false`
// is an ordinary unitary operation.
//
// The spatial part is composed, not re-declared: there is one definition of a
// rotation-plus-translation, and `spatial` is it.
struct MagneticSymmetryOperation {
  SymmetryOperation spatial{};
  bool time_reversal = false;
};

} // namespace cppcrystal

#pragma GCC visibility pop
