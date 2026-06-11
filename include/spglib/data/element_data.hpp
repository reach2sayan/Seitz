#pragma once

#include <spglib/data/element_tables.hpp>

#include <numbers>
#include <optional>
#include <string_view>

// Accessors over the curated element table (generated element_tables.hpp:
// covalent radii from Cordero et al. 2008). Atom "types" in this codebase are
// atomic numbers, so these are looked up directly by type. The covalent radius
// is the primary datum; the atomic volume used to size a generated cell is
// derived from it as a sphere volume (the PyXtal radius-based estimate), so no
// second table is sourced.
namespace spglib::data {

// True iff `z` is a tabulated atomic number (1..kNumElements).
[[nodiscard]] constexpr bool is_known_element(int z) noexcept {
  return z >= 1 && z <= kNumElements;
}

// Single-bond covalent radius (angstrom) of atomic number `z`; std::nullopt for
// an untabulated element (no magic sentinel).
[[nodiscard]] constexpr std::optional<double>
covalent_radius(int z) noexcept {
  if (!is_known_element(z)) {
    return std::nullopt;
  }
  return kCovalentRadii[static_cast<std::size_t>(z - 1)];
}

// Estimated atomic volume (cubic angstrom) of atomic number `z`, the volume of
// a sphere of its covalent radius. std::nullopt for an untabulated element.
[[nodiscard]] constexpr std::optional<double>
atomic_volume(int z) noexcept {
  auto const r = covalent_radius(z);
  if (!r) {
    return std::nullopt;
  }
  return 4.0 / 3.0 * std::numbers::pi * (*r) * (*r) * (*r);
}

// Chemical symbol of atomic number `z`; std::nullopt for an untabulated element.
[[nodiscard]] constexpr std::optional<std::string_view>
element_symbol(int z) noexcept {
  if (!is_known_element(z)) {
    return std::nullopt;
  }
  return kElementSymbols[static_cast<std::size_t>(z - 1)];
}

// Atomic number of a chemical symbol (case-sensitive, e.g. "Na" -> 11);
// std::nullopt if the symbol is not tabulated. Backed by a Boost.Bimap built
// once from the generated table (see src/data/element_data.cpp).
[[nodiscard]] std::optional<int> atomic_number(std::string_view symbol);

} // namespace spglib::data
