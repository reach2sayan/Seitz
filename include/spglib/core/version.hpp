#pragma once

namespace spglib {

struct Version {
  int major;
  int minor;
  int patch;
};

// This port's own version.
inline constexpr Version kVersion{0, 1, 0};

// The reference spglib release whose behaviour this port targets and validates
// against (see the oracle test setup).
inline constexpr Version kReferenceSpglibVersion{2, 7, 0};

// "major.minor.patch" of this port.
[[nodiscard]] char const *version_string() noexcept;

} // namespace spglib
