#include <cppcrystal/data/element_data.hpp>

#include <boost/bimap.hpp>

#include <string>

namespace cppcrystal::data {

namespace {

// Symbol <-> atomic number, built once from the generated table. A Boost.Bimap
// rather than two parallel maps so the two directions can never drift apart.
// (Only the symbol -> number direction is currently exposed; the reverse is
// already covered by the constexpr element_symbol().)
using SymbolMap = boost::bimap<std::string, int>;

[[nodiscard]] SymbolMap const &symbol_map() {
  static SymbolMap const map = [] {
    SymbolMap m;
    for (int z = 1; z <= kNumElements; ++z) {
      m.insert(SymbolMap::value_type(
          std::string(kElementSymbols[static_cast<std::size_t>(z - 1)]), z));
    }
    return m;
  }();
  return map;
}

} // namespace

std::optional<int> atomic_number(std::string_view symbol) {
  auto const &map = symbol_map();
  auto const it = map.left.find(std::string(symbol));
  if (it == map.left.end()) {
    return std::nullopt;
  }
  return it->second;
}

} // namespace cppcrystal::data
