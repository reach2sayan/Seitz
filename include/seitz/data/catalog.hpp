#pragma once

#include <seitz/core/keys.hpp>

#include <array>
#include <cstddef>

// The compile-time catalogs of the built-in databases: one class over a family
// policy, keyed by a validated key (HallNumber, UniNumber).
#pragma GCC visibility push(default)

namespace seitz::data {

// A family policy supplies
// (a) the decoded row type,
// (b) how many there are,
// (c) the key that addresses them,
// (d) the decode from the generated raw table.
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

// One catalog per family, decoded from the generated tables at compile time.
// `decode()` yields the row array, not a Catalog: a policy naming
// Catalog<Self> in its own body would make CatalogFamily<Self> self-dependent.
template <CatalogFamily Family>
inline constexpr Catalog<Family> kCatalog{Family::decode()};

} // namespace seitz::data

#pragma GCC visibility pop
