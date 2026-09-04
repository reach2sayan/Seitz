#pragma once

// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal {

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

} // namespace cppcrystal

#pragma GCC visibility pop
