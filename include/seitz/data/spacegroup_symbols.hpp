#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/data/catalog.hpp>
#include <seitz/data/detail/lookup.hpp>
#include <seitz/data/spg_database.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string_view>

// The symbol -> Hall setting direction of the database, which nothing else
// needed until a CIF arrived carrying `P 21/c` instead of a symop loop. Its own
// header rather than part of spg_database.hpp: three translation units want it
// and it costs a few thousand consteval normalisations.
//
// Symbols in the wild are written every way at once -- `P 21/c`, `P 1 21/c 1`,
// `P2_1/c`, `P2(1)/c` are the same group -- so both sides of the comparison are
// normalised to a canonical key: spaces, underscores and parentheses dropped,
// and everything after a ':' taken as the setting the writer asked for
// (`R -3 m :H`). Case is kept, because it is what separates a space group's `P`
// from a layer group's `p`.

#pragma GCC visibility push(default)

namespace seitz::data {

namespace detail {

// Longest normalised symbol in either table is 13 characters (the longest Hall
// symbol, 11); the setting suffix is at most a few.
inline constexpr std::size_t kSymbolCapacity = 16;

// How many keys one catalog row contributes: the short symbol, the full
// symbol, up to three ` = ` alternatives of the international symbol, and for
// the monoclinic rows the full symbol with its lone `1` axes removed
// (`P 1 2_1/n 1` -> `P21/n`). Absent kinds normalise to the empty key, which
// no query can equal.
inline constexpr int kSymbolKinds = 6;
inline constexpr std::size_t kSymbolBuckets = 4096;

// A normalised symbol, as a value a consteval table can hold.
struct SymbolKey {
  std::array<char, kSymbolCapacity> text{};
  std::uint8_t size = 0;

  [[nodiscard]] constexpr std::string_view view() const noexcept {
    return {text.data(), size};
  }
  [[nodiscard]] constexpr bool empty() const noexcept { return size == 0; }
  [[nodiscard]] friend constexpr bool operator==(SymbolKey const &a,
                                                 SymbolKey const &b) noexcept {
    return a.view() == b.view();
  }
};

struct NormalizedSymbol {
  SymbolKey key;
  SymbolKey setting;
};

[[nodiscard]] constexpr bool is_dropped(char c) noexcept {
  return c == ' ' || c == '_' || c == '(' || c == ')';
}

// Append `c`, or report the overflow. A symbol longer than the capacity is one
// no row carries, so its caller abandons the whole key rather than truncating
// it into something that could collide with a real one.
[[nodiscard]] constexpr bool push(SymbolKey &key, char c) noexcept {
  if (key.size >= kSymbolCapacity) {
    return false;
  }
  key.text[key.size++] = c;
  return true;
}

// `symbol` split into its canonical key and the setting the writer asked for.
// `drop_lone_ones` removes whole `1` tokens, which is what turns a monoclinic
// full symbol into the short one everybody writes.
[[nodiscard]] constexpr NormalizedSymbol
normalize_symbol(std::string_view symbol, bool drop_lone_ones = false) noexcept {
  NormalizedSymbol out{};
  if (auto const colon = symbol.find(':'); colon != std::string_view::npos) {
    for (char const c : symbol.substr(colon + 1)) {
      if (!is_dropped(c) && !push(out.setting, c)) {
        return {};
      }
    }
    symbol = symbol.substr(0, colon);
  }
  while (!symbol.empty()) {
    auto const space = symbol.find(' ');
    std::string_view const token = symbol.substr(0, space);
    symbol = space == std::string_view::npos ? std::string_view{}
                                             : symbol.substr(space + 1);
    if (token.empty() || (drop_lone_ones && token == "1")) {
      continue;
    }
    for (char const c : token) {
      if (!is_dropped(c) && !push(out.key, c)) {
        return {};
      }
    }
  }
  return out;
}

// The `n`-th ` = ` alternative of an international symbol, empty past the last.
[[nodiscard]] constexpr std::string_view alternative(std::string_view symbol,
                                                     int n) noexcept {
  constexpr std::string_view separator = " = ";
  for (int i = 0; i < n; ++i) {
    auto const at = symbol.find(separator);
    if (at == std::string_view::npos) {
      return {};
    }
    symbol = symbol.substr(at + separator.size());
  }
  auto const at = symbol.find(separator);
  return at == std::string_view::npos ? symbol : symbol.substr(0, at);
}

// The pre-1983 cubic spelling modernised: `Pm3m` is `Pm-3m`, and CIF files
// written before the change are still in every database. Empty when the symbol
// carries no such `3`, so the retry is skipped.
[[nodiscard]] constexpr SymbolKey modernized_cubic(SymbolKey const &key) noexcept {
  auto const text = key.view();
  SymbolKey out{};
  bool changed = false;
  for (std::size_t i = 0; i < text.size(); ++i) {
    if (text[i] == '3' && i > 0 &&
        std::string_view{"mnabcde"}.find(text[i - 1]) !=
            std::string_view::npos) {
      if (!push(out, '-')) {
        return {};
      }
      changed = true;
    }
    if (!push(out, text[i])) {
      return {};
    }
  }
  return changed ? out : SymbolKey{};
}

[[nodiscard]] constexpr int bucket_of(SymbolKey const &key) noexcept {
  std::uint32_t hash = 2166136261U;
  for (char const c : key.view()) {
    hash ^= static_cast<std::uint8_t>(c);
    hash *= 16777619U;
  }
  return static_cast<int>(hash % kSymbolBuckets);
}

// An id of the Hermann-Mauguin index addresses one (row, kind) pair.
[[nodiscard]] constexpr int hm_row_of(int id) noexcept {
  return (id - 1) / kSymbolKinds + 1;
}
[[nodiscard]] constexpr int hm_kind_of(int id) noexcept {
  return (id - 1) % kSymbolKinds;
}

template <GroupFamily F>
[[nodiscard]] constexpr SpacegroupType const &row_at(int row) noexcept {
  return kCatalog<SpacegroupFamily<F>>.rows[static_cast<std::size_t>(row) - 1];
}

template <GroupFamily F>
[[nodiscard]] constexpr SymbolKey hm_key_of(int id) noexcept {
  SpacegroupType const &type = row_at<F>(hm_row_of(id));
  switch (hm_kind_of(id)) {
  case 0:
    return normalize_symbol(type.international_short).key;
  case 1:
    return normalize_symbol(type.international_full).key;
  case 2:
  case 3:
  case 4:
    return normalize_symbol(alternative(type.international, hm_kind_of(id) - 2))
        .key;
  default:
    // Monoclinic only: `P312` and `P321` are different trigonal groups, so
    // dropping lone `1`s there would merge them.
    return type.number >= 3 && type.number <= 15
               ? normalize_symbol(type.international_full, true).key
               : SymbolKey{};
  }
}

template <GroupFamily F>
inline constexpr auto kHmSymbolIndex =
    bucket_index<static_cast<std::size_t>(hall_settings(F)) * kSymbolKinds,
                 kSymbolBuckets>([](int id) { return bucket_of(hm_key_of<F>(id)); });

// A Hall symbol is already one per setting, so its index is keyed by the row
// itself rather than by a (row, kind) pair.
template <GroupFamily F>
[[nodiscard]] constexpr SymbolKey hall_key_of(int row) noexcept {
  return normalize_symbol(row_at<F>(row).hall_symbol).key;
}

template <GroupFamily F>
inline constexpr auto kHallSymbolIndex =
    bucket_index<static_cast<std::size_t>(hall_settings(F)), kSymbolBuckets>(
        [](int row) { return bucket_of(hall_key_of<F>(row)); });

// An unset setting takes whatever the first matching row offers; a set one has
// to prefix the row's choice, so `:H` and `:R` pick different rows of the same
// group and `:b` picks the first of `b1`, `b2`, `b3`.
template <GroupFamily F>
[[nodiscard]] constexpr bool choice_accepts(int row,
                                            SymbolKey const &setting) noexcept {
  return setting.empty() || row_at<F>(row).choice.starts_with(setting.view());
}

// The first (lowest-index, therefore default) setting whose `key_of` equals
// `key` and whose choice accepts `setting`. Bucket entries are ascending in id,
// so the scan meets the rows in catalog order.
template <GroupFamily F, class Index, class KeyOf, class RowOf>
[[nodiscard]] constexpr std::optional<HallNumber>
find_setting(Index const &index, KeyOf key_of, RowOf row_of,
             SymbolKey const &key, SymbolKey const &setting) noexcept {
  if (key.empty()) {
    return std::nullopt;
  }
  for (int const id : index[bucket_of(key)]) {
    if (key_of(id) == key && choice_accepts<F>(row_of(id), setting)) {
      return HallNumber::of(F, row_of(id));
    }
  }
  return std::nullopt;
}

} // namespace detail

// The Hall setting a Hermann-Mauguin symbol names, in whichever of the
// spellings the tables carry: short, full, an alternative of the international
// symbol, or a monoclinic full symbol with its lone axes dropped. A trailing
// `:choice` selects among a group's settings; without one the default (lowest
// Hall number) wins. std::nullopt when nothing matches.
//
// Layer-group short symbols are not unique -- `p2` is both layer group 3 and 8
// -- and there the lower group number wins.
template <GroupFamily F>
[[nodiscard]] constexpr std::optional<HallNumber>
hall_from_hm_symbol(std::string_view symbol) noexcept {
  auto const query = detail::normalize_symbol(symbol);
  auto const search = [&](detail::SymbolKey const &key) {
    return detail::find_setting<F>(detail::kHmSymbolIndex<F>,
                                   detail::hm_key_of<F>, detail::hm_row_of, key,
                                   query.setting);
  };
  if (auto const found = search(query.key)) {
    return found;
  }
  return search(detail::modernized_cubic(query.key));
}

// The Hall setting a Hall symbol names (`-P 2ybc`, `-F 4 2 3`). Hall symbols
// are already one per setting, so no choice suffix is needed -- one is still
// honoured if written.
template <GroupFamily F>
[[nodiscard]] constexpr std::optional<HallNumber>
hall_from_hall_symbol(std::string_view symbol) noexcept {
  auto const query = detail::normalize_symbol(symbol);
  return detail::find_setting<F>(detail::kHallSymbolIndex<F>,
                                 detail::hall_key_of<F>,
                                 [](int row) { return row; }, query.key,
                                 query.setting);
}

// ---- compile-time guards ---------------------------------------------------

namespace detail {
// The four spellings of one monoclinic group all reach the same setting.
constexpr bool monoclinic_spellings_agree() {
  auto const a = hall_from_hm_symbol<GroupFamily::space>("P 21/c");
  return a.has_value() &&
         a == hall_from_hm_symbol<GroupFamily::space>("P 1 21/c 1") &&
         a == hall_from_hm_symbol<GroupFamily::space>("P2_1/c") &&
         a == hall_from_hm_symbol<GroupFamily::space>("P2(1)/c") &&
         spacegroup_type(*a).number == 14;
}
} // namespace detail

static_assert(detail::monoclinic_spellings_agree());
// The pre-1983 cubic spelling and the modern one are the same group.
static_assert(
    spacegroup_type(*hall_from_hm_symbol<GroupFamily::space>("Pm3m")).number ==
    221);
static_assert(hall_from_hm_symbol<GroupFamily::space>("Pm3m") ==
              hall_from_hm_symbol<GroupFamily::space>("Pm-3m"));
// The rhombohedral and hexagonal settings of R-3m are one group, two rows.
static_assert(hall_from_hm_symbol<GroupFamily::space>("R -3 m :R") !=
              hall_from_hm_symbol<GroupFamily::space>("R -3 m :H"));
static_assert(
    spacegroup_type(*hall_from_hm_symbol<GroupFamily::space>("R -3 m :R"))
        .number == 166);
static_assert(
    spacegroup_type(*hall_from_hm_symbol<GroupFamily::space>("R -3 m :H"))
        .number == 166);
// A Hall symbol names its setting directly.
static_assert(
    spacegroup_type(*hall_from_hall_symbol<GroupFamily::space>("-F 4 2 3"))
        .number == 225);
// An alternative orthorhombic setting keeps its choice.
static_assert(
    spacegroup_type(*hall_from_hm_symbol<GroupFamily::space>("Pbnm")).choice ==
    "cab");
static_assert(hall_from_hm_symbol<GroupFamily::layer>("p 2/m 1 1").has_value());
static_assert(!hall_from_hm_symbol<GroupFamily::space>("Xyz").has_value());
static_assert(!hall_from_hm_symbol<GroupFamily::space>("").has_value());
static_assert(!hall_from_hall_symbol<GroupFamily::space>("nonsense").has_value());

} // namespace seitz::data

#pragma GCC visibility pop
