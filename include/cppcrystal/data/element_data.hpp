#pragma once

#include <cppcrystal/data/element_tables.hpp>

#include <numbers>
#include <optional>
#include <string_view>

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

[[nodiscard]] std::optional<int> atomic_number(std::string_view symbol);

} // namespace cppcrystal::data
