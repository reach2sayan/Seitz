#include <seitz/core/symmetry_operation.hpp>

#include "math/integer_matrix.hpp"

#include <boost/algorithm/string/join.hpp>
#include <boost/parser/parser.hpp>

#include <array>
#include <cmath>
#include <cstdlib>
#include <format>
#include <optional>
#include <ranges>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

// The Jones-faithful coordinate triplet, both directions. The grammar is
// private to this TU: Boost::parser is PRIVATE and no installed header names it.
namespace seitz {

namespace {

namespace bp = boost::parser;

// One summand: an axis multiple or a bare number. The sign is its own field
// because the grammar reads it before knowing which of the two follows.
struct RotTerm {
  std::optional<char> sign;
  std::optional<unsigned> coefficient;
  int axis;
};
struct TransTerm {
  std::optional<char> sign;
  double value;
};
using Term = std::variant<RotTerm, TransTerm>;
using Coordinate = std::vector<Term>;
struct Triplet {
  Coordinate a;
  Coordinate b;
  Coordinate c;
};

bp::symbols<int> const axis_symbol{
    {"x", 0}, {"y", 1}, {"z", 2}, {"X", 0}, {"Y", 1}, {"Z", 2}};
auto const sign = bp::char_("+-");

// Fraction first, so "1/2" is one term; backtracks to a decimal on no '/'.
bp::rule<struct number_tag, double> const number = "number";
auto const number_def =
    bp::transform([](std::tuple<unsigned, unsigned> const &pq) {
      auto const [p, q] = pq;
      return q == 0 ? 0.0 : static_cast<double>(p) / static_cast<double>(q);
    })[bp::uint_ >> '/' >> bp::uint_] |
    bp::double_;

bp::rule<struct rot_term_tag, RotTerm> const rot_term = "rotation term";
auto const rot_term_def = -sign >> -bp::uint_ >> axis_symbol;

bp::rule<struct trans_term_tag, TransTerm> const trans_term =
    "translation term";
auto const trans_term_def = -sign >> number;

// Rotation first: on "1/2+x" it takes the 1, fails at '/', and backtracks.
bp::rule<struct coordinate_tag, Coordinate> const coordinate = "coordinate";
auto const coordinate_def = +(rot_term | trans_term);

// separate[]: without it the three coordinates' vectors merge into one.
bp::rule<struct triplet_tag, Triplet> const triplet = "xyz triplet";
auto const triplet_def =
    bp::separate[coordinate >> ',' >> coordinate >> ',' >> coordinate];

// The macro's generated functions take an unused `dont_assign`; it expands
// here, so Clang blames this file despite parser.hpp being a SYSTEM include.
#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
BOOST_PARSER_DEFINE_RULES(number, rot_term, trans_term, coordinate, triplet);
#if defined(__clang__)
#pragma clang diagnostic pop
#endif

[[nodiscard]] constexpr int sign_of(std::optional<char> c) noexcept {
  return c == '-' ? -1 : 1;
}

constexpr std::string_view kAxisNames = "xyz";

// The denominators a translation can have exactly (Hall settings carry 3/8).
constexpr std::array kDenominators = {1, 2, 3, 4, 6, 8, 12};
constexpr double kFractionPrec = 1e-6;

// `t` as "+p/q" when exact, else a short decimal. `signed_`: does anything
// precede it, so is a leading '+' needed.
void append_translation(std::string &out, double t, bool signed_) {
  if (std::abs(t) < kFractionPrec) {
    return;
  }
  auto const exact = [t](int q) {
    return std::abs(t - std::round(t * q) / q) <= kFractionPrec;
  };
  if (auto const q = std::ranges::find_if(kDenominators, exact);
      q != kDenominators.end()) {
    auto const p = static_cast<long>(std::lround(t * *q));
    if (p >= 0 && signed_) {
      out += '+';
    }
    out += *q == 1 ? std::format("{}", p) : std::format("{}/{}", p, *q);
    return;
  }
  if (t >= 0.0 && signed_) {
    out += '+';
  }
  out += std::format("{:.6g}", t);
}

// The axis terms, then the translation. "0" for an empty row, so it reads back.
[[nodiscard]] std::string coordinate_text(SymmetryOperation const &op,
                                          Index row) {
  std::string out;
  for (auto const [axis, coefficient] :
       std::views::iota(Index{0}, Index{3}) |
           std::views::transform([&](Index a) {
             return std::pair{a, op.rotation(row, a)};
           }) |
           std::views::filter([](auto const &t) { return t.second != 0; })) {
    if (coefficient < 0) {
      out += '-';
    } else if (!out.empty()) {
      out += '+';
    }
    if (std::abs(coefficient) != 1) {
      out += std::format("{}", std::abs(coefficient));
    }
    out += kAxisNames[static_cast<std::size_t>(axis)];
  }
  append_translation(out, op.translation[row], !out.empty());
  return out.empty() ? std::string{"0"} : out;
}

} // namespace

std::string to_xyz(SymmetryOperation const &op) {
  auto const coordinates =
      std::views::iota(Index{0}, Index{3}) |
      std::views::transform([&](Index row) { return coordinate_text(op, row); });
  return boost::algorithm::join(
      std::vector<std::string>{std::from_range, coordinates}, ",");
}

Result<SymmetryOperation> from_xyz(std::string_view text) {
  auto const invalid = [&] {
    return leaf::new_error(e_invalid_xyz{std::string{text}});
  };

  Triplet parsed;
  if (!bp::parse(text, triplet, bp::ws, parsed)) {
    return invalid();
  }

  SymmetryOperation op{.rotation = Matrix3i::Zero(),
                       .translation = Vector3d::Zero()};
  for (auto const [row, coordinate_terms] :
       std::views::enumerate(std::array{&parsed.a, &parsed.b, &parsed.c})) {
    for (Term const &term : *coordinate_terms) {
      std::visit(
          [&op, r = static_cast<Index>(row)](auto const &t) {
            if constexpr (std::same_as<std::remove_cvref_t<decltype(t)>,
                                       RotTerm>) {
              op.rotation(r, t.axis) +=
                  sign_of(t.sign) * static_cast<int>(t.coefficient.value_or(1));
            } else {
              op.translation[r] += sign_of(t.sign) * t.value;
            }
          },
          term);
    }
  }

  if (std::abs(math::determinant(math::as_rows(op.rotation))) != 1) {
    return invalid();
  }
  return op;
}

} // namespace seitz
