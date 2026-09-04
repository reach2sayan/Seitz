#pragma once

#include <cppcrystal/core/keys.hpp>

#include <array>
#include <cstddef>

// The compile-time catalogs of the built-in databases. One class parameterised
// by a family policy replaces the three hand-written catalog structs, and the
// index-0 sentinel row they each carried is gone: a Catalog is keyed by a
// validated key (HallNumber, UniNumber), so every lookup is in range by
// construction and there is nothing to fall back to.
// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal::data {

// A family policy supplies the decoded row type, how many there are, the key
// that addresses them, and the decode from the generated raw table.
template <class Family>
concept CatalogFamily = requires {
  typename Family::Row;
  typename Family::Key;
  { Family::count } -> std::convertible_to<std::size_t>;
};

// Rows addressed by a validated 1-based key. `rows[i]` holds key index i + 1.
template <CatalogFamily Family> struct Catalog {
  using Row = typename Family::Row;
  using Key = typename Family::Key;
  static constexpr std::size_t size = Family::count;

  std::array<Row, size> rows{};

  // Total: Key cannot name an out-of-range entry.
  [[nodiscard]] constexpr Row const &operator[](Key key) const noexcept {
    return rows[static_cast<std::size_t>(Family::index_of(key)) - 1];
  }

  [[nodiscard]] constexpr auto begin() const noexcept { return rows.begin(); }
  [[nodiscard]] constexpr auto end() const noexcept { return rows.end(); }
};

// The one catalog per family, decoded from the generated tables at compile
// time. `decode()` yields the row array rather than a Catalog: a family policy
// that named Catalog<Self> in its own body would make CatalogFamily<Self>
// depend on itself.
template <CatalogFamily Family>
inline constexpr Catalog<Family> kCatalog{Family::decode()};

} // namespace cppcrystal::data

#pragma GCC visibility pop
