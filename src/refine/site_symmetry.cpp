#include "refine/site_symmetry.hpp"

#include "core/position_index.hpp"
#include "data/sitesym_database.hpp"
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include <seitz/core/operation_set.hpp>
#include <seitz/core/periodicity.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

// Site symmetry (3D space-group path). Determines the exact (symmetrized)
// location of each conventional-cell atom from its site symmetry (Grosse-
// Kunstleve & Adams 2002), groups equivalent atoms, and looks up the Wyckoff
// letter + site-symmetry symbol from the database.
namespace seitz::refine {

namespace {

constexpr int kNumAttempt = 5;
constexpr double kIncreaseRate = 1.05;

// Average the site-symmetry operations that fix `position` to pin it onto its
// exact special position.
[[nodiscard]] Vector3d set_exact_location(Vector3d const &position,
                                          Operations const &conv_sym,
                                          Matrix3d const &lattice,
                                          double symprec,
                                          CellPeriodicity const &periodicity) {
  Matrix3d sum_rot = Matrix3d::Zero();
  Vector3d sum_trans = Vector3d::Zero();
  int num = 0;
  for (auto const &op : conv_sym) {
    Vector3d const pos = op.apply(position);
    if (coincident(pos, position, lattice, symprec, periodicity)) {
      sum_rot += op.rotation.cast<double>();
      Vector3d const wrap = math::round_to_int(pos - position).cast<double>();
      sum_trans += op.translation - wrap;
      ++num;
    }
  }
  sum_rot /= num;
  sum_trans /= num;
  return sum_rot * position + sum_trans;
}

struct Equivalent {
  Vector3d position;
  int representative;
};

// Place every atom either at the exact special position (independent atoms)
// or as the image of a representative. Once an independent atom is placed, its
// orbit under the conventional operations claims every later atom it lands on,
// in (representative, operation) order; the first claim wins. Wyckoff letters
// and symbols are filled in later by set_wyckoff_labels.
[[nodiscard]] ExactPositions get_exact_positions(Cell const &conv_prim,
                                                 Operations const &conv_sym,
                                                 double symprec) {
  CellPeriodicity const &periodicity = conv_prim.periodicity();
  PositionIndex const index(conv_prim, symprec);
  PositionIndex::Scratch scratch; // reused by every query in the loop below
  auto const n = static_cast<std::size_t>(conv_prim.size());
  std::vector<std::optional<Equivalent>> claimed(n);
  ExactPositions out;
  out.reserve(n);

  for (int i = 0; i < conv_prim.size(); ++i) {
    if (auto const &image = claimed[static_cast<std::size_t>(i)]) {
      out.push_back({.position = image->position,
                     .equivalent_atom = image->representative});
      continue;
    }
    Vector3d const exact =
        set_exact_location(conv_prim.position(i), conv_sym,
                           conv_prim.lattice().matrix(), symprec, periodicity);
    out.push_back({.position = exact, .equivalent_atom = i});
    for (auto const &op : conv_sym) {
      Vector3d const mapped = op.apply(exact);
      for (int const k : index.matches(mapped, conv_prim.type(i), scratch)) {
        auto &slot = claimed[static_cast<std::size_t>(k)];
        if (k > i && !slot) {
          slot = Equivalent{math::wrap_to_unit_cell(mapped), i};
        }
      }
    }
  }
  return out;
}

struct WyckoffLabel {
  int letter;
  // A view into the static site-symmetry table, never a copy.
  std::string_view symbol;
};

// Match `position`'s orbit against the database Wyckoff positions of the Hall
// setting; the first consistent one gives the letter and site-symmetry symbol.
// A candidate is consistent when, for some orbit point, the number of orbit
// points coinciding with it AND fixed by the candidate's site-symmetry
// generator, times the candidate multiplicity, equals the group order.
[[nodiscard]] std::optional<WyckoffLabel>
get_wyckoff_notation(Vector3d const &position, Operations const &conv_sym,
                     int ref_multiplicity, Matrix3d const &lattice,
                     HallNumber hall, double symprec,
                     CellPeriodicity const &periodicity) {
  std::vector<Vector3d> const orbit(
      std::from_range, conv_sym | std::views::transform([&](auto const &op) {
                         return op.apply(position);
                       }));

  Positions const orbit_positions = to_positions(orbit);
  Types const orbit_types(orbit.size(), 0);
  PositionIndex const index(orbit_positions, orbit_types, lattice, symprec,
                            periodicity);
  PositionIndex::Scratch scratch;

  auto const group_order = static_cast<int>(conv_sym.size());
  data::WyckoffRange const range = data::wyckoff_indices(hall);
  // A candidate is consistent when some orbit point j has exactly
  // group_order / multiplicity fixed points among the orbit points coinciding
  // with it. Counted from the fixed side: coincidence is symmetric, so every
  // fixed point k adds one to the count of each j coinciding with k -- which
  // needs one index query per fixed point instead of one per orbit point, and
  // none at all for the usual wrong candidate, whose generator fixes nothing.
  // Fixed points at identical coordinates share their query.
  std::vector<int> fixed;
  fixed.reserve(orbit.size());
  std::vector<int> counts(orbit.size());
  for (int wi = 0; wi < range.count; ++wi) {
    data::WyckoffCoordinate const wc =
        data::wyckoff_coordinate(range.start + wi);
    if (wc.multiplicity != ref_multiplicity) {
      continue;
    }
    auto const fixes = [&](int k) {
      Vector3d const &point = orbit[static_cast<std::size_t>(k)];
      Vector3d const mapped =
          wc.rotation.cast<double>() * point + wc.translation;
      return coincident(point, mapped, lattice, symprec, periodicity);
    };
    fixed.clear();
    std::ranges::copy(std::views::iota(0, static_cast<int>(orbit.size())) |
                          std::views::filter(fixes),
                      std::back_inserter(fixed));
    if (fixed.empty()) {
      continue;
    }
    auto const coordinates = [&](int k) {
      auto const &p = orbit[static_cast<std::size_t>(k)];
      return std::array{p[0], p[1], p[2]};
    };
    std::ranges::sort(fixed, {}, coordinates);
    std::ranges::fill(counts, 0);
    for (auto const chunk : fixed | std::views::chunk_by([&](int a, int b) {
                              return coordinates(a) == coordinates(b);
                            })) {
      auto const weight = static_cast<int>(std::ranges::size(chunk));
      for (int const j : index.matches(
               orbit[static_cast<std::size_t>(chunk.front())], scratch)) {
        counts[static_cast<std::size_t>(j)] += weight;
      }
    }
    bool const consistent = std::ranges::any_of(counts, [&](int count) {
      return count * wc.multiplicity == group_order;
    });
    if (consistent) {
      // The database stores Wyckoff positions in reverse order (g f e ... a).
      return WyckoffLabel{range.count - wi - 1,
                          data::site_symmetry_symbol(range.start + wi)};
    }
  }
  return std::nullopt;
}

// Assign a Wyckoff letter + symbol to each independent atom, then propagate to
// the equivalent atoms. False if any lookup failed.
[[nodiscard]] bool set_wyckoff_labels(ExactPositions &atoms,
                                      Cell const &conv_prim,
                                      Operations const &conv_sym,
                                      int num_pure_trans, HallNumber hall,
                                      double symprec) {
  CellPeriodicity const &periodicity = conv_prim.periodicity();

  std::vector<int> nums_equiv(atoms.size(), 0);
  for (auto const &atom : atoms) {
    ++nums_equiv[static_cast<std::size_t>(atom.equivalent_atom)];
  }

  for (auto [i, atom] : atoms | std::views::enumerate) {
    if (atom.equivalent_atom != i) {
      continue;
    }
    auto const label = get_wyckoff_notation(
        atom.position, conv_sym,
        nums_equiv[static_cast<std::size_t>(i)] * num_pure_trans,
        conv_prim.lattice().matrix(), hall, symprec, periodicity);
    if (!label) {
      return false;
    }
    atom.wyckoff = label->letter;
    atom.site_symmetry_symbol = label->symbol;
  }

  for (auto [i, atom] : atoms | std::views::enumerate) {
    if (atom.equivalent_atom != i) {
      auto const &rep = atoms[static_cast<std::size_t>(atom.equivalent_atom)];
      atom.wyckoff = rep.wyckoff;
      atom.site_symmetry_symbol = rep.site_symmetry_symbol;
    }
  }
  return true;
}

} // namespace

std::optional<ExactPositions> exact_positions(Cell const &conv_prim,
                                              Operations const &conv_sym,
                                              int num_pure_trans,
                                              HallNumber hall, double symprec) {
  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    ExactPositions exact = get_exact_positions(conv_prim, conv_sym, tolerance);
    if (set_wyckoff_labels(exact, conv_prim, conv_sym, num_pure_trans, hall,
                           symprec)) {
      return exact;
    }
    tolerance *= kIncreaseRate;
  }
  return std::nullopt;
}

} // namespace seitz::refine
