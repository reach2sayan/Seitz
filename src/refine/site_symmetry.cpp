#include <cppcrystal/refine/site_symmetry.hpp>

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/position_index.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <algorithm>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

// Site symmetry (3D space-group path). Determines the exact (symmetrized)
// location of each conventional-cell atom from its site symmetry (Grosse-
// Kunstleve & Adams 2002), groups equivalent atoms, and looks up the Wyckoff
// letter + site-symmetry symbol from the database.
namespace cppcrystal::refine {

namespace {

constexpr int kNumAttempt = 5;
constexpr double kIncreaseRate = 1.05;

// Average the site-symmetry operations that fix `position` to pin it onto its
// exact special position.
[[nodiscard]] Vector3d set_exact_location(Vector3d const &position,
                                          SymmetryOperations const &conv_sym,
                                          Matrix3d const &lattice,
                                          double symprec,
                                          std::optional<int> aperiodic_axis) {
  Matrix3d sum_rot = Matrix3d::Zero();
  Vector3d sum_trans = Vector3d::Zero();
  int num = 0;
  for (auto const &op : conv_sym) {
    Vector3d const pos = op.apply(position);
    if (is_overlap(pos, position, lattice, symprec, aperiodic_axis)) {
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
[[nodiscard]] ExactPositions
get_exact_positions(Cell const &conv_prim, SymmetryOperations const &conv_sym,
                    double symprec) {
  std::optional<int> const aperiodic_axis = conv_prim.aperiodic_axis();
  PositionIndex const index(conv_prim, symprec);
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
                           conv_prim.lattice(), symprec, aperiodic_axis);
    out.push_back({.position = exact, .equivalent_atom = i});
    for (auto const &op : conv_sym) {
      Vector3d const mapped = op.apply(exact);
      for (int const k : index.matches(mapped, conv_prim.type(i))) {
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
  std::string symbol;
};

// Match `position`'s orbit against the database Wyckoff positions of the Hall
// setting; the first consistent one gives the letter and site-symmetry symbol.
// A candidate is consistent when, for some orbit point, the number of orbit
// points coinciding with it AND fixed by the candidate's site-symmetry
// generator, times the candidate multiplicity, equals the group order.
[[nodiscard]] std::optional<WyckoffLabel>
get_wyckoff_notation(Vector3d const &position, SymmetryOperations const &conv_sym,
                     int ref_multiplicity, Matrix3d const &lattice,
                     int hall_number, double symprec,
                     std::optional<int> aperiodic_axis) {
  std::vector<Vector3d> orbit;
  orbit.reserve(conv_sym.size());
  std::ranges::transform(conv_sym, std::back_inserter(orbit),
                         [&](auto const &op) { return op.apply(position); });

  // Coincidence classes of the orbit, found once through an index over it.
  CellPeriodicity const periodicity =
      periodicity_from_aperiodic_axis(aperiodic_axis);
  Positions const orbit_positions = to_positions(orbit);
  Types const orbit_types(orbit.size(), 0);
  PositionIndex const index(BucketGeometry::of(lattice, symprec, periodicity),
                            orbit_positions, orbit_types, lattice, symprec,
                            periodicity);
  std::vector<std::vector<int>> coincident;
  coincident.reserve(orbit.size());
  for (auto const &point : orbit) {
    coincident.push_back(std::ranges::to<std::vector<int>>(index.matches(point)));
  }

  auto const group_order = static_cast<int>(conv_sym.size());
  data::WyckoffRange const range = data::wyckoff_indices(hall_number);
  for (int wi = 0; wi < range.count; ++wi) {
    data::WyckoffCoordinate const wc =
        data::wyckoff_coordinate(range.start + wi);
    if (wc.multiplicity != ref_multiplicity) {
      continue;
    }
    // Which orbit points the candidate's site-symmetry generator fixes.
    std::vector<bool> fixed(orbit.size());
    for (auto const [k, point] : orbit | std::views::enumerate) {
      Vector3d const mapped = wc.rotation.cast<double>() * point + wc.translation;
      fixed[static_cast<std::size_t>(k)] =
          is_overlap(point, mapped, lattice, symprec, aperiodic_axis);
    }
    bool const consistent =
        std::ranges::any_of(coincident, [&](std::vector<int> const &cls) {
          auto const fixed_in_class = std::ranges::count_if(
              cls, [&](int k) { return fixed[static_cast<std::size_t>(k)]; });
          return static_cast<int>(fixed_in_class) * wc.multiplicity ==
                 group_order;
        });
    if (consistent) {
      // The database stores Wyckoff positions in reverse order (g f e ... a).
      return WyckoffLabel{
          range.count - wi - 1,
          std::string(data::site_symmetry_symbol(range.start + wi))};
    }
  }
  return std::nullopt;
}

// Assign a Wyckoff letter + symbol to each independent atom, then propagate to
// the equivalent atoms. False if any lookup failed.
[[nodiscard]] bool set_wyckoff_labels(ExactPositions &atoms,
                                      Cell const &conv_prim,
                                      SymmetryOperations const &conv_sym,
                                      int num_pure_trans, int hall_number,
                                      double symprec) {
  std::optional<int> const aperiodic_axis = conv_prim.aperiodic_axis();

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
        conv_prim.lattice(), hall_number, symprec, aperiodic_axis);
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

std::optional<ExactPositions>
exact_positions(Cell const &conv_prim, SymmetryOperations const &conv_sym,
                int num_pure_trans, int hall_number, double symprec) {
  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    ExactPositions exact = get_exact_positions(conv_prim, conv_sym, tolerance);
    if (set_wyckoff_labels(exact, conv_prim, conv_sym, num_pure_trans,
                           hall_number, symprec)) {
      return exact;
    }
    tolerance *= kIncreaseRate;
  }
  return std::nullopt;
}

} // namespace cppcrystal::refine
