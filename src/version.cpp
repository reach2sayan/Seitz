#include <seitz/core/version.hpp>

namespace seitz {

// SEITZ_VERSION_STRING comes from the configured version.hpp, which CMake
// writes from project(Seitz VERSION ...). Out of line so that the string has
// one address across the library rather than one per translation unit.
char const *version_string() noexcept { return SEITZ_VERSION_STRING; }

} // namespace seitz
