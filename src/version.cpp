#include <seitz/core/version.hpp>

namespace seitz {

// SEITZ_VERSION_STRING comes from the configured version.hpp. Out of line so
// the string has one address across the library, not one per TU.
char const *version_string() noexcept { return SEITZ_VERSION_STRING; }

} // namespace seitz
