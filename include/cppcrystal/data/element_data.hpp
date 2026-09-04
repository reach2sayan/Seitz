#pragma once

#include <cppcrystal/data/element_tables.hpp>

#include <algorithm>
#include <array>
#include <numbers>
#include <optional>
#include <string_view>
#include <utility>

// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal::data {

// True iff `z` is a tabulated atomic number (1..kNumElements).
[[nodiscard]] constexpr bool is_known_element(int z) noexcept {
  return z >= 1 && z <= kNumElements;
}

// Single-bond covalent radius (angstrom) of atomic number `z`; std::nullopt for
// an untabulated element (no magic sentinel).
[[nodiscard]] constexpr std::optional<double> covalent_radius(int z) noexcept {
  if (!is_known_element(z)) {
    return std::nullopt;
  }
  return kCovalentRadii[static_cast<std::size_t>(z - 1)];
}

[[nodiscard]] constexpr std::optional<double> atomic_volume(int z) noexcept {
  auto const r = covalent_radius(z);
  if (!r) {
    return std::nullopt;
  }
  return 4.0 / 3.0 * std::numbers::pi * (*r) * (*r) * (*r);
}

// Chemical symbol of atomic number `z`; std::nullopt for an untabulated
// element.
[[nodiscard]] constexpr std::optional<std::string_view>
element_symbol(int z) noexcept {
  if (!is_known_element(z)) {
    return std::nullopt;
  }
  return kElementSymbols[static_cast<std::size_t>(z - 1)];
}

namespace detail {
// (symbol, atomic number) sorted by symbol, for a binary search at compile
// time or run time.
using SymbolEntry = std::pair<std::string_view, int>;
inline constexpr auto kSymbolsByName = [] {
  std::array<SymbolEntry, kNumElements> t{};
  for (int z = 1; z <= kNumElements; ++z) {
    t[static_cast<std::size_t>(z - 1)] = {
        kElementSymbols[static_cast<std::size_t>(z - 1)], z};
  }
  std::ranges::sort(t, {}, &SymbolEntry::first);
  return t;
}();
} // namespace detail

// Atomic number of a chemical symbol (case-sensitive); std::nullopt for an
// unknown symbol.
[[nodiscard]] constexpr std::optional<int>
atomic_number(std::string_view symbol) noexcept {
  auto const it = std::ranges::lower_bound(detail::kSymbolsByName, symbol, {},
                                           &detail::SymbolEntry::first);
  if (it == detail::kSymbolsByName.end() || it->first != symbol) {
    return std::nullopt;
  }
  return it->second;
}

static_assert(atomic_number("H") == 1);
static_assert(atomic_number(*element_symbol(kNumElements)) == kNumElements);
static_assert(!atomic_number("Xx"));

} // namespace cppcrystal::data

#pragma GCC visibility pop
