#include <seitz/io/cif.hpp>

#include <seitz/core/fractional.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/data/element_data.hpp>
#include <seitz/data/spacegroup_symbols.hpp>
#include <seitz/data/spg_database.hpp>

#include "core/position_index.hpp"
#include "math/lattice_parameters.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/join.hpp>
// Declared ahead of the parser headers: Boost.Parser's tracing prints every
// sub-parser through a qualified call, so a user parser's overload has to be
// visible where those templates are defined, not merely where they are used.
#include <iosfwd>
namespace seitz::io {
template <class Token> struct TokenParser;
}
namespace boost::parser::detail {
template <class Context, class Token>
void print_parser(Context const &, seitz::io::TokenParser<Token> const &,
                  std::ostream &os, int components = 0);
}

#include <boost/parser/parser.hpp>
#include <boost/range/join.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <locale>
#include <cmath>
#include <cstdint>
#include <format>
#include <iterator>
#include <numbers>
#include <optional>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

// Grammar and raw aggregates are private to this TU: Boost::parser is a
// PRIVATE dependency and no installed header names it.
namespace seitz::io {

namespace {

namespace bp = boost::parser;

using Iterator = std::string_view::const_iterator;

struct Item {
  std::string tag;
  std::string value;
};
struct Loop {
  Iterator where{}; // the `loop_` keyword, so a ragged loop reports its line
  std::vector<std::string> tags;
  std::vector<std::string> values; // row-major
};
using Entry = std::variant<Loop, Item>;
struct Block {
  std::string name;
  std::vector<Entry> entries;
};

// ---- token parsers ----------------------------------------------------------
//
// Boost.Parser charges per parser invocation, not per character: every
// sub-parser (even inside lexeme[]) enters the skip machinery and sets up a
// parse context before it looks at a byte, so `*(char_ - ws)` over a value
// cost ~30 ns a character and the parse was two thirds of a CIF read, and
// `!no_case[keyword]` case-folded every value through Unicode tables. Each
// CIF token is therefore one leaf parser, built the way the library's own
// char_ and string_ leaves are (the two `call` overloads, detail::append for
// the attribute, a print_parser hook for tracing) so it composes with >>, |,
// lexeme[] and rules like any other. A Token policy says how a token opens,
// where it ends, how it closes, and whether the text is acceptable; the
// grammar below is still the combinator one.
//
// Whitespace and line ends are the ASCII ones bp::ws and bp::eol match on
// `char` input: space, tab through carriage return; \n, \r (\r\n as one),
// \v, \f.
[[nodiscard]] constexpr bool is_blank(char c) noexcept {
  return c == ' ' || ('\t' <= c && c <= '\r');
}
[[nodiscard]] constexpr bool is_eol(char c) noexcept {
  return c == '\n' || c == '\r' || c == '\v' || c == '\f';
}
[[nodiscard]] constexpr char ascii_lower(char c) noexcept {
  return 'A' <= c && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
}

// Reserved words a value may not begin with (case-insensitively); they end
// whatever list they follow. `save_` here means a dictionary's save frames
// stop the parse -- out of scope.
[[nodiscard]] bool starts_with_keyword(Iterator first, Iterator last) {
  std::string_view const text(first, last);
  return std::ranges::any_of(
      std::array{std::string_view{"data_"}, std::string_view{"loop_"},
                 std::string_view{"save_"}, std::string_view{"global_"},
                 std::string_view{"stop_"}},
      [&](std::string_view keyword) {
        return text.size() >= keyword.size() &&
               std::ranges::equal(keyword, text.substr(0, keyword.size()), {},
                                  {}, ascii_lower);
      });
}

} // namespace

template <class Token> struct TokenParser {
  template <class Iter, class Sentinel, class Context, class SkipParser>
  std::string call(Iter &first, Sentinel last, Context const &context,
                   SkipParser const &skip, bp::detail::flags flags,
                   bool &success) const {
    std::string out;
    call(first, last, context, skip, flags, success, out);
    return out;
  }
  template <class Iter, class Sentinel, class Context, class SkipParser,
            class Attribute>
  void call(Iter &first, Sentinel last, Context const &context,
            SkipParser const &, bp::detail::flags flags, bool &success,
            Attribute &out) const {
    Iter it = first;
    if (!Token::open(it, last, context)) {
      success = false;
      return;
    }
    Iter const begin = it;
    while (it != last && !Token::at_end(it, last)) {
      ++it;
    }
    Iter const end = it;
    if (!Token::accept(begin, end) || !Token::close(it, last)) {
      success = false;
      return;
    }
    // Appended as a range, not through detail::append, whose range form
    // push_backs one character at a time.
    if constexpr (!bp::detail::is_nope_v<Attribute>) {
      if (bp::detail::gen_attrs(flags)) {
        auto const at = out.size();
        out.insert(out.end(), begin, end);
        if constexpr (Token::lowercase) {
          std::ranges::transform(out.begin() + static_cast<std::ptrdiff_t>(at),
                                 out.end(),
                                 out.begin() + static_cast<std::ptrdiff_t>(at),
                                 ascii_lower);
        }
      }
    }
    first = it;
  }
};

namespace {

// The default policy: no delimiters, any text, ends at whitespace.
struct Bare {
  static constexpr bool lowercase = false;
  static bool open(auto &, auto, auto const &) { return true; }
  static bool at_end(auto it, auto) { return is_blank(*it); }
  static bool accept(auto, auto) { return true; }
  static bool close(auto &, auto) { return true; }
};
// The rest of the line, for comments.
struct LineRest : Bare {
  static bool at_end(auto it, auto) { return is_eol(*it); }
};
// `_tag`, underscore included, lowercased: a CIF tag is ASCII, and a Turkish
// one dotless-Is it.
struct Tag : Bare {
  static constexpr bool lowercase = true;
  static bool open(auto &it, auto last, auto const &) {
    return it != last && *it == '_';
  }
  static bool accept(auto begin, auto end) { return end - begin > 1; }
};
// An unquoted value: may not start like a tag, a quote, a comment, a text
// field or a reserved character, and may not be a keyword.
struct Unquoted : Bare {
  static bool open(auto &it, auto last, auto const &) {
    return it != last && !is_blank(*it) &&
           std::string_view{"_'\"#;$[]"}.find(*it) == std::string_view::npos;
  }
  static bool accept(auto begin, auto end) {
    return !starts_with_keyword(begin, end);
  }
};
// A quoted value; the closing quote is the one followed by whitespace or the
// end, so 'it's here' holds.
template <char Quote> struct Quoted : Bare {
  static bool open(auto &it, auto last, auto const &) {
    if (it == last || *it != Quote) {
      return false;
    }
    ++it;
    return true;
  }
  static bool at_end(auto it, auto last) {
    return *it == Quote && (std::next(it) == last || is_blank(*std::next(it)));
  }
  static bool close(auto &it, auto last) {
    if (it == last) {
      return false;
    }
    ++it;
    return true;
  }
};
// A text field: ';' opening a line through the next line end followed by ';'.
struct TextField : Bare {
  static bool open(auto &it, auto last, auto const &context) {
    bool const at_column_1 =
        it == bp::_begin(context) || *std::prev(it) == '\n';
    if (!at_column_1 || it == last || *it != ';') {
      return false;
    }
    ++it;
    return true;
  }
  static bool at_end(auto it, auto last) {
    if (!is_eol(*it)) {
      return false;
    }
    auto next = std::next(it);
    if (*it == '\r' && next != last && *next == '\n') {
      ++next;
    }
    return next != last && *next == ';';
  }
  static bool close(auto &it, auto last) {
    if (it == last) {
      return false;
    }
    if (*it == '\r' && std::next(it) != last && *std::next(it) == '\n') {
      ++it;
    }
    ++it; // the line end
    ++it; // the ';'
    return true;
  }
};

} // namespace
} // namespace seitz::io

// The tracing hook declared above the parser headers; nothing here traces,
// but the sequence parser instantiates it regardless.
namespace boost::parser::detail {
template <class Context, class Token>
void print_parser(Context const &, seitz::io::TokenParser<Token> const &,
                  std::ostream &os, int) {
  os << "token";
}
} // namespace boost::parser::detail

namespace seitz::io {
namespace {

bp::parser_interface<TokenParser<Bare>> const bare;
bp::parser_interface<TokenParser<LineRest>> const line_rest;
bp::parser_interface<TokenParser<Tag>> const cif_tag;
bp::parser_interface<TokenParser<Unquoted>> const unquoted;
bp::parser_interface<TokenParser<Quoted<'\''>>> const squoted;
bp::parser_interface<TokenParser<Quoted<'"'>>> const dquoted;
bp::parser_interface<TokenParser<TextField>> const text_field;

auto const comment = bp::lit('#') >> bp::omit[line_rest];
auto const skipper = bp::ws | comment;

// Skipped once, before the token; the alternatives are single leaves.
auto const value = bp::lexeme[text_field | squoted | dquoted | unquoted];

// separate[]: two adjacent string leaves would otherwise merge into one.
bp::rule<struct item_tag, Item> const item = "tag/value item";
auto const item_def = bp::separate[cif_tag >> value];

// separate[]: without it the tag and value lists merge into one vector.
bp::rule<struct loop_tag, Loop> const loop = "loop";
auto const loop_def = bp::separate[bp::transform([](auto const &r) {
                                     return r.begin();
                                   })[bp::raw[bp::no_case[bp::lit("loop_")]]] >>
                                   +cif_tag >> +value];

bp::rule<struct block_tag, Block> const data_block = "data block";
auto const data_block_def =
    bp::lexeme[bp::no_case[bp::lit("data_")] >> bare] >> *(loop | item);

// The trailing `*skipper` is what lets the caller check `first == last`:
// blocks alone stop at the last value, leaving the file's final newline.
bp::rule<struct document_tag, std::vector<Block>> const document = "CIF";
auto const document_def = *data_block >> bp::omit[*skipper];

// The macro's generated functions take an unused `dont_assign`; it expands
// here, so Clang blames this file despite parser.hpp being a SYSTEM include.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
BOOST_PARSER_DEFINE_RULES(item, loop, data_block, document);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

[[nodiscard]] leaf::error_id syntax_error(Iterator begin, Iterator at) {
  auto const position = bp::find_line_position(begin, at);
  return leaf::new_error(
      e_cif_syntax{.line = position.line_number + 1,
                   .column = position.column_number + 1});
}

} // namespace

std::optional<std::span<std::string const>>
CifBlock::column(std::string_view tag) const {
  auto const it = columns.find(tag);
  if (it == columns.end()) {
    return std::nullopt;
  }
  return std::span<std::string const>{it->second};
}

std::optional<std::string_view> CifBlock::value(std::string_view tag) const {
  auto const rows = column(tag);
  if (!rows || rows->empty()) {
    return std::nullopt;
  }
  std::string_view const first = rows->front();
  if (first == "?" || first == ".") {
    return std::nullopt;
  }
  return first;
}

namespace {

// ---- the value grammar the semantic layer needs -----------------------------

// esd stripped: `1.234(5)` is 1.234, and an unclosed `0.25(` is still 0.25.
[[nodiscard]] std::optional<double> number_of(std::string_view text) {
  double out = 0.0;
  auto const number =
      bp::lexeme[bp::double_ >>
                 -bp::omit['(' >> *bp::digit >> -bp::lit(')')]];
  if (!bp::parse(text, number, bp::ws, out)) {
    return std::nullopt;
  }
  return out;
}

// The element a label or type symbol names: leading letters, two first so
// `Cl1` is chlorine not carbon, then one for `O1W`, `Fe3+`, `Na+`.
[[nodiscard]] std::optional<int> species_of(std::string_view symbol) {
  std::string letters{
      std::from_range,
      symbol | std::views::take_while([](unsigned char c) {
        return std::isalpha(c) != 0;
      })};
  // The hydrogen isotopes neutron structures name are no element of their own.
  if (letters == "D" || letters == "T") {
    return data::atomic_number("H");
  }
  auto const titled = [&](std::size_t n) {
    std::string s =
        boost::algorithm::to_lower_copy(letters.substr(0, n),
                                        std::locale::classic());
    s.front() = static_cast<char>(
        std::toupper(static_cast<unsigned char>(s.front())));
    return s;
  };
  for (std::size_t n : {std::size_t{2}, std::size_t{1}}) {
    if (letters.size() >= n) {
      if (auto const z = data::atomic_number(titled(n))) {
        return z;
      }
    }
  }
  return std::nullopt;
}

// ---- the asymmetric unit ----------------------------------------------------

struct CifSite {
  std::string label;
  int type = 0;
  Vector3d position{Vector3d::Zero()};
  double occupancy = 1.0;
  int representative = 0; // lowest site index coincident with this one
};

[[nodiscard]] leaf::error_id missing(std::string_view tag) {
  return leaf::new_error(e_cif_missing{std::string{tag}});
}

[[nodiscard]] Result<double> scalar_of(CifBlock const &block,
                                       std::string_view tag) {
  auto const text = block.value(tag);
  if (!text) {
    return missing(tag);
  }
  auto const number = number_of(*text);
  if (!number) {
    return leaf::new_error(e_cif_missing{std::string{tag}},
                           e_message{"not a number: " + std::string{*text}});
  }
  return *number;
}

// An optional column, padded to `rows` so every site column zips.
[[nodiscard]] std::vector<std::string>
column_or_blank(CifBlock const &block, std::string_view tag,
                std::size_t rows) {
  auto const column = block.column(tag);
  if (!column || column->size() != rows) {
    return std::vector<std::string>(rows);
  }
  return {column->begin(), column->end()};
}

[[nodiscard]] Result<Lattice> lattice_of(CifBlock const &block) {
  constexpr double kDegree = std::numbers::pi / 180.0;
  BOOST_LEAF_AUTO(a, scalar_of(block, "_cell_length_a"));
  BOOST_LEAF_AUTO(b, scalar_of(block, "_cell_length_b"));
  BOOST_LEAF_AUTO(c, scalar_of(block, "_cell_length_c"));
  BOOST_LEAF_AUTO(alpha, scalar_of(block, "_cell_angle_alpha"));
  BOOST_LEAF_AUTO(beta, scalar_of(block, "_cell_angle_beta"));
  BOOST_LEAF_AUTO(gamma, scalar_of(block, "_cell_angle_gamma"));

  Matrix3d const basis = math::lattice_from_parameters(
      a, b, c, alpha * kDegree, beta * kDegree, gamma * kDegree);
  auto lattice = Lattice::from_basis(basis);
  if (!lattice) {
    return leaf::new_error(e_invalid_lattice{basis.determinant()});
  }
  return std::move(*lattice);
}

[[nodiscard]] Result<std::vector<CifSite>> sites_of(CifBlock const &block) {
  auto const x = block.column("_atom_site_fract_x");
  auto const y = block.column("_atom_site_fract_y");
  auto const z = block.column("_atom_site_fract_z");
  if (!x) {
    return missing("_atom_site_fract_x");
  }
  if (!y) {
    return missing("_atom_site_fract_y");
  }
  if (!z) {
    return missing("_atom_site_fract_z");
  }
  if (y->size() != x->size() || z->size() != x->size()) {
    return leaf::new_error(e_cif_missing{"_atom_site_fract_y"},
                           e_message{"the coordinate columns differ in length"});
  }

  auto const rows = x->size();
  auto const labels = column_or_blank(block, "_atom_site_label", rows);
  auto const symbols = column_or_blank(block, "_atom_site_type_symbol", rows);
  auto const occupancies =
      column_or_blank(block, "_atom_site_occupancy", rows);

  std::vector<CifSite> sites;
  sites.reserve(rows);
  for (auto const [row, columns] : std::views::enumerate(
           std::views::zip(*x, *y, *z, labels, symbols, occupancies))) {
    auto const &[sx, sy, sz, label, symbol, occupancy] = columns;
    auto const coordinate = std::array{number_of(sx), number_of(sy),
                                       number_of(sz)};
    if (!std::ranges::all_of(coordinate, &std::optional<double>::has_value)) {
      return leaf::new_error(
          e_cif_missing{"_atom_site_fract_x"},
          e_message{"a coordinate is not a number: " + sx + ' ' + sy + ' ' +
                    sz});
    }

    CifSite site;
    site.label = label.empty() ? std::format("site{}", row) : label;
    site.position = Vector3d{*coordinate[0], *coordinate[1], *coordinate[2]};
    // "Not stated" is a full site when reading permissively.
    site.occupancy = number_of(occupancy).value_or(1.0);

    // '?' or '.' says the file has none, not that the element is called "?".
    bool const stated =
        !symbol.empty() && symbol != "?" && symbol != ".";
    std::string_view const named = stated ? symbol : site.label;
    auto const type = species_of(named);
    if (!type) {
      return leaf::new_error(e_unknown_element{std::string{named}});
    }
    site.type = *type;
    sites.push_back(std::move(site));
  }
  return sites;
}

// Sites sharing a position collapsed to the majority species. Fills in each
// site's representative; returns the survivors, ascending, and the records.
[[nodiscard]] std::pair<std::vector<int>, std::vector<OccupancyCollapse>>
collapse_shared_sites(std::vector<CifSite> &sites, Lattice const &lattice,
                      double symprec) {
  // Non-owning index, so both arrays are named locals. One type throughout:
  // sharing a position is geometric, and a mixed site is where species differ.
  std::vector<Vector3d> const points{
      std::from_range, sites | std::views::transform(&CifSite::position)};
  Positions const positions = to_positions(points);
  Types const one_type(sites.size(), 0);
  PositionIndex const index{positions, one_type, lattice.matrix(), symprec,
                            all_periodic()};

  for (auto const [row, site] : std::views::enumerate(sites)) {
    auto matches = index.matches(site.position);
    auto const first = std::ranges::begin(matches);
    site.representative =
        first == std::ranges::end(matches) ? static_cast<int>(row) : *first;
  }

  std::vector<int> order{
      std::from_range,
      std::views::iota(0, static_cast<int>(sites.size()))};
  std::ranges::stable_sort(order, {}, [&](int i) {
    return sites[static_cast<std::size_t>(i)].representative;
  });

  std::vector<int> kept;
  std::vector<OccupancyCollapse> collapsed;
  auto const occupancy_of = [&](int i) {
    return sites[static_cast<std::size_t>(i)].occupancy;
  };
  auto const label_of = [&](int i) {
    return sites[static_cast<std::size_t>(i)].label;
  };
  for (auto chunk : order | std::views::chunk_by([&](int a, int b) {
                      return sites[static_cast<std::size_t>(a)].representative ==
                             sites[static_cast<std::size_t>(b)].representative;
                    })) {
    int const winner = std::ranges::max(chunk, {}, occupancy_of);
    kept.push_back(winner);
    if (std::ranges::distance(chunk) > 1 || occupancy_of(winner) < 1.0) {
      // Named first: MSVC cannot parse a capturing lambda inside a designated
      // initializer's nested braces, and loses the enclosing scope after it.
      std::vector<std::string> dropped(
          std::from_range, chunk | std::views::filter([&](int i) {
                             return i != winner;
                           }) | std::views::transform(label_of));
      collapsed.push_back(OccupancyCollapse{.kept = label_of(winner),
                                            .occupancy = occupancy_of(winner),
                                            .dropped = std::move(dropped)});
    }
  }
  std::ranges::sort(kept);
  return {std::move(kept), std::move(collapsed)};
}

// ---- the symmetry the block states ------------------------------------------

constexpr auto kSymopTags = std::array<std::string_view, 2>{
    "_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"};
constexpr auto kHallTags = std::array<std::string_view, 2>{
    "_space_group_name_hall", "_symmetry_space_group_name_hall"};
constexpr auto kHmTags = std::array<std::string_view, 2>{
    "_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m"};
constexpr auto kNumberTags = std::array<std::string_view, 2>{
    "_space_group_it_number", "_symmetry_int_tables_number"};

[[nodiscard]] std::optional<std::string_view>
first_stated(CifBlock const &block, std::span<std::string_view const> tags) {
  for (std::string_view const tag : tags) {
    if (auto const stated = block.value(tag)) {
      return stated;
    }
  }
  return std::nullopt;
}

[[nodiscard]] std::optional<HallNumber>
hall_of_number(std::string_view text) {
  auto const number = number_of(text);
  return number ? data::default_hall<GroupFamily::space>(
                      static_cast<int>(std::lround(*number)))
                : std::nullopt;
}

// The three ways a block can name its setting, most specific first.
struct NamedSource {
  std::span<std::string_view const> tags;
  std::optional<HallNumber> (*resolve)(std::string_view);
};
constexpr auto kNamedSources = std::array{
    NamedSource{kHallTags, &data::hall_from_hall_symbol<GroupFamily::space>},
    NamedSource{kHmTags, &data::hall_from_hm_symbol<GroupFamily::space>},
    NamedSource{kNumberTags, &hall_of_number}};

// The setting the block names and the tables know; nullopt otherwise. For
// the loop path, where the name only confirms the listed operations, so an
// unknown symbol (`Aba2`, a non-standard setting) is no error.
[[nodiscard]] std::optional<HallNumber>
resolved_setting(CifBlock const &block) {
  for (auto const &[tags, resolve] : kNamedSources) {
    auto const stated = first_stated(block, tags);
    if (auto const hall = stated ? resolve(*stated) : std::nullopt) {
      return hall;
    }
  }
  return std::nullopt;
}

// The setting the block names outright; nullopt when it names none, an error
// when it names one the tables cannot resolve.
[[nodiscard]] Result<std::optional<HallNumber>>
named_setting(CifBlock const &block) {
  for (auto const &[tags, resolve] : kNamedSources) {
    auto const stated = first_stated(block, tags);
    if (!stated) {
      continue;
    }
    auto const hall = resolve(*stated);
    if (!hall) {
      return leaf::new_error(e_unknown_spacegroup_symbol{std::string{*stated}});
    }
    return hall;
  }
  return std::optional<HallNumber>{};
}

// Whether `listed` is exactly the operation set of `hall`: the same count,
// and every listed operation present with the same rotation and a coincident
// translation (modulo lattice translations, within symprec). n^2 over at
// most 192 operations; the search this spares walks 530 candidate settings.
[[nodiscard]] bool lists_setting(std::span<SymmetryOperation const> listed,
                                 HallNumber hall, Lattice const &lattice,
                                 double symprec) {
  Operations const &database = data::operations_from_database(hall);
  if (listed.size() != database.size()) {
    return false;
  }
  return std::ranges::all_of(listed, [&](SymmetryOperation const &op) {
    return std::ranges::any_of(database, [&](SymmetryOperation const &ref) {
      return ref.rotation == op.rotation &&
             coincident(ref.translation, op.translation, lattice.matrix(),
                        symprec, all_periodic());
    });
  });
}

// The operations to expand with, and their setting. A symop loop wins -- it is
// what the file asserts; failing to recover its Hall setting is not fatal.
// A setting the block also names is taken when the loop lists exactly its
// operations, which is what a well-formed file does; otherwise the setting is
// searched for from the operations alone.
[[nodiscard]] Result<std::pair<std::vector<SymmetryOperation>,
                               std::optional<HallNumber>>>
symmetry_of(CifBlock const &block, Lattice const &lattice,
            Tolerance const &tol) {
  for (std::string_view const tag : kSymopTags) {
    auto const column = block.column(tag);
    if (!column) {
      continue;
    }
    std::vector<SymmetryOperation> operations;
    operations.reserve(column->size());
    for (std::string const &text : *column) {
      BOOST_LEAF_AUTO(operation, from_xyz(text));
      operations.push_back(operation);
    }
    std::optional<HallNumber> const named = resolved_setting(block);
    std::optional<HallNumber> hall =
        named && lists_setting(operations, *named, lattice, tol.symprec)
            ? named
            : std::nullopt;
    if (!hall) {
      if (auto const match =
              Operations{operations}.spacegroup(lattice.matrix(), tol)) {
        hall = match->hall;
      }
    }
    return std::pair{std::move(operations), hall};
  }

  BOOST_LEAF_AUTO(hall, named_setting(block));
  if (!hall) {
    return std::pair{std::vector<SymmetryOperation>{SymmetryOperation{}},
                     std::optional<HallNumber>{}};
  }
  Operations const &database = data::operations_from_database(*hall);
  return std::pair{
      std::vector<SymmetryOperation>{database.begin(), database.end()}, hall};
}

} // namespace

Result<std::vector<CifBlock>> parse_cif(std::string_view text) {
  std::vector<Block> raw;
  auto first = text.begin();
  auto const last = text.end();
  bool const matched = bp::prefix_parse(first, last, document, skipper, raw);
  if (!matched || first != last) {
    return syntax_error(text.begin(), first);
  }

  std::vector<CifBlock> blocks;
  blocks.reserve(raw.size());
  for (Block &parsed : raw) {
    CifBlock out;
    out.name = std::move(parsed.name);
    for (Entry &entry : parsed.entries) {
      if (Item *scalar = std::get_if<Item>(&entry); scalar != nullptr) {
        out.columns.insert_or_assign(
            std::move(scalar->tag),
            std::vector<std::string>{std::move(scalar->value)});
        continue;
      }
      Loop &table = std::get<Loop>(entry);
      auto const width = table.tags.size();
      if (table.values.size() % width != 0) {
        return syntax_error(text.begin(), table.where);
      }
      for (auto const [offset, name] : std::views::enumerate(table.tags)) {
        out.columns.insert_or_assign(
            name, std::vector<std::string>{
                      std::from_range,
                      table.values |
                          std::views::drop(static_cast<std::size_t>(offset)) |
                          std::views::stride(width)});
      }
      out.loops.push_back(std::move(table.tags));
    }
    blocks.push_back(std::move(out));
  }
  return blocks;
}


Result<CifStructure> structure_of(CifBlock const &block, Tolerance tol) {
  BOOST_LEAF_AUTO(lattice, lattice_of(block));
  BOOST_LEAF_AUTO(sites, sites_of(block));
  auto [kept, collapsed] = collapse_shared_sites(sites, lattice, tol.symprec);
  BOOST_LEAF_AUTO(symmetry, symmetry_of(block, lattice, tol));
  auto const &[operations, hall] = symmetry;

  // Every image of every kept site, in (site, operation) order; an image is
  // placed unless an earlier placed one coincides with it -- against
  // everything placed so far, not just this site's images: files do list a
  // row and its own image, and collapse_shared_sites only sees rows already
  // coincident before the symmetry is applied. One index over all the
  // candidates answers "an earlier placed coincident image?" in order, which
  // is the same rule as placing one by one against the placed list.
  Positions images(static_cast<Index>(kept.size() * operations.size()), 3);
  std::vector<int> site_of;
  site_of.reserve(static_cast<std::size_t>(images.rows()));
  for (int const index : kept) {
    CifSite const &site = sites[static_cast<std::size_t>(index)];
    for (SymmetryOperation const &operation : operations) {
      images.row(static_cast<Index>(site_of.size())) =
          math::wrap_to_unit_cell(operation.apply(site.position)).transpose();
      site_of.push_back(index);
    }
  }
  Types const untyped(site_of.size(), 0); // coincidence ignores the species
  PositionIndex const index(images, untyped, lattice.matrix(), tol.symprec,
                            all_periodic());
  PositionIndex::Scratch scratch;
  std::vector<std::uint8_t> placed(site_of.size(), 0);
  std::vector<Vector3d> positions;
  Types types;
  std::vector<std::string> labels;
  for (auto const [i, site_index] : std::views::enumerate(site_of)) {
    auto const earlier = index.first_match(
        images.row(i).transpose(), 0, scratch,
        [&](int j) { return j < i && placed[static_cast<std::size_t>(j)]; });
    if (earlier) {
      continue;
    }
    placed[static_cast<std::size_t>(i)] = 1;
    CifSite const &site = sites[static_cast<std::size_t>(site_index)];
    positions.push_back(images.row(i).transpose());
    types.push_back(site.type);
    labels.push_back(site.label);
  }
  if (positions.empty()) {
    return leaf::new_error(e_empty_cell{});
  }

  return CifStructure{.name = block.name,
                      .cell = Cell{lattice, to_positions(positions),
                                   std::move(types)},
                      .hall = hall,
                      .labels = std::move(labels),
                      .collapsed = std::move(collapsed)};
}

Result<std::vector<CifStructure>> read_cif(std::string_view text,
                                           Tolerance tol) {
  BOOST_LEAF_AUTO(blocks, parse_cif(text));
  std::vector<CifStructure> structures;
  for (CifBlock const &block : blocks) {
    if (!block.column("_cell_length_a")) {
      continue;
    }
    BOOST_LEAF_AUTO(structure, structure_of(block, tol));
    structures.push_back(std::move(structure));
  }
  return structures;
}


namespace {

// ---- the writer -------------------------------------------------------------

// One `_atom_site` row: an orbit's representative, named and counted.
struct Row {
  std::string label;
  std::string_view symbol;
  int multiplicity = 1;
  char wyckoff = 'a';
  Vector3d position{Vector3d::Zero()};
};

// Quoted, since every symbol worth writing carries spaces.
[[nodiscard]] std::string quoted_value(std::string_view text) {
  return std::format("'{}'", text);
}

// `Fe1`, `Fe2`, ... : one counter per species, in the order the rows come.
[[nodiscard]] std::vector<Row>
labelled(std::vector<Row> rows) {
  boost::container::flat_map<std::string_view, int> counter;
  for (Row &row : rows) {
    row.label = std::format("{}{}", row.symbol, ++counter[row.symbol]);
  }
  return rows;
}

[[nodiscard]] std::string_view symbol_of(int type) {
  // `X` is what CIF readers treat as an unknown scatterer.
  return data::element_symbol(type).value_or("X");
}

[[nodiscard]] std::string render(std::string_view name, Lattice const &lattice,
                                 HallNumber hall,
                                 std::span<SymmetryOperation const> operations,
                                 std::vector<Row> const &rows) {
  data::SpacegroupType const &type = data::spacegroup_type(hall);
  constexpr double kDegree = 180.0 / std::numbers::pi;
  Matrix3d const metric = lattice.metric();

  // The underscores of `P 1 2_1/c 1` are a table convention, not the symbol.
  std::string hm{std::from_range,
                 type.international_full |
                     std::views::filter([](char c) { return c != '_'; })};
  // A non-default setting must say so, or a symbol-only reader lands on the
  // group's default (`R -3 2/m` alone is hexagonal, never rhombohedral).
  if (data::default_hall<GroupFamily::space>(type.number) != hall) {
    hm += std::format(" :{}", type.choice);
  }

  std::string out = std::format("data_{}\n", name);
  out += std::format("_space_group_name_H-M_alt      {}\n", quoted_value(hm));
  out += std::format("_space_group_name_Hall         {}\n",
                     quoted_value(type.hall_symbol));
  out += std::format("_space_group_IT_number         {}\n", type.number);
  // One sequence, so the row format is written once. Each angle is named by
  // the axis it is opposite: alpha = b^c, beta = c^a, gamma = a^b.
  constexpr auto kCellTags =
      std::array{"_cell_length_a",    "_cell_length_b",   "_cell_length_c",
                 "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma"};
  auto const angle = [&](int i) {
    return math::metric_angle(metric, (i + 1) % 3, (i + 2) % 3) * kDegree;
  };
  std::array const lengths{math::metric_length(metric, 0),
                           math::metric_length(metric, 1),
                           math::metric_length(metric, 2)};
  std::array const angles{angle(0), angle(1), angle(2)};
  for (auto const [tag, parameter] :
       std::views::zip(kCellTags, boost::range::join(lengths, angles))) {
    out += std::format("{:<30} {:.6f}\n", tag, parameter);
  }

  out += "\nloop_\n_space_group_symop_id\n_space_group_symop_operation_xyz\n";
  for (auto const [i, operation] : std::views::enumerate(operations)) {
    out += std::format("{} {}\n", i + 1, quoted_value(to_xyz(operation)));
  }

  out += "\nloop_\n_atom_site_label\n_atom_site_type_symbol\n"
         "_atom_site_symmetry_multiplicity\n_atom_site_Wyckoff_symbol\n"
         "_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n"
         "_atom_site_occupancy\n";
  for (Row const &row : rows) {
    out += std::format("{} {} {} {} {:.6f} {:.6f} {:.6f} {:.4f}\n", row.label,
                       row.symbol, row.multiplicity, row.wyckoff,
                       row.position[0], row.position[1], row.position[2], 1.0);
  }
  return out;
}

} // namespace

std::string write_cif(Cell const &cell, std::string_view name) {
  auto const rows = labelled(std::vector<Row>{
      std::from_range,
      cell.atoms() | std::views::transform([](auto const &atom) {
        auto const &[position, type] = atom;
        return Row{.label = {},
                   .symbol = symbol_of(type),
                   .multiplicity = 1,
                   .wyckoff = 'a',
                   .position = position};
      })});
  constexpr auto kP1 = *HallNumber::of(GroupFamily::space, 1);
  return render(name, cell.lattice(), kP1,
                std::array{SymmetryOperation{}}, rows);
}

Result<std::string> write_cif(analysis::SymmetryAnalyzer const &analyzer,
                              std::string_view name) {
  BOOST_LEAF_AUTO(hall, analyzer.hall());
  BOOST_LEAF_AUTO(standardized, analyzer.standardized_cell());

  // Re-analysed with the setting fixed: the Wyckoff assignment and operations
  // must be the ones valid in the basis being written, not the input cell's.
  auto const restated = analysis::SymmetryAnalyzer::from_cell(
      standardized, analyzer.tolerance(), hall);
  BOOST_LEAF_AUTO(sites, restated.sites());
  BOOST_LEAF_AUTO(operations, restated.operations());

  auto const rows = labelled(std::vector<Row>{
      std::from_range,
      std::views::enumerate(sites) |
          std::views::filter([](auto const &indexed) {
            auto const &[index, site] = indexed;
            return static_cast<int>(index) == site.equivalent_atom;
          }) |
          std::views::transform([&](auto const &indexed) {
            auto const &[index, site] = indexed;
            return Row{
                .label = {},
                .symbol = symbol_of(standardized.type(static_cast<Index>(index))),
                .multiplicity = static_cast<int>(std::ranges::count(
                    sites, static_cast<int>(index),
                    &analysis::Site::equivalent_atom)),
                .wyckoff = static_cast<char>('a' + site.wyckoff),
                .position = standardized.position(static_cast<Index>(index))};
          })});

  return render(name, standardized.lattice(), hall, operations.span(), rows);
}

} // namespace seitz::io
