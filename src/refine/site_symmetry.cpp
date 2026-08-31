#include <cppcrystal/refine/site_symmetry.hpp>

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <algorithm>
#include <iterator>
#include <optional>
#include <ranges>

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

// If atom `i` is the image of an already-placed independent atom under a
// conventional operation, return that placement.
[[nodiscard]] std::optional<Equivalent>
find_equivalent(int i, std::vector<int> const &indep_atoms,
                ExactPositions const &placed, Cell const &conv_prim,
                SymmetryOperations const &conv_sym, double symprec) {
  std::optional<int> const aperiodic_axis = conv_prim.aperiodic_axis();
  for (int const j : indep_atoms) {
    for (auto const &op : conv_sym) {
      Vector3d const pos =
          op.apply(placed[static_cast<std::size_t>(j)].position);
      if (is_overlap_same_type(pos, conv_prim.position(i), conv_prim.type(j),
                               conv_prim.type(i), conv_prim.lattice(), symprec,
                               aperiodic_axis)) {
        return Equivalent{math::wrap_to_unit_cell(pos), j};
      }
    }
  }
  return std::nullopt;
}

// Place every atom either at the exact special position (independent atoms) or
// as the image of a representative (equivalent atoms). Wyckoff letters and
// symbols are filled in later by set_wyckoff_labels.
[[nodiscard]] ExactPositions
get_exact_positions(Cell const &conv_prim, SymmetryOperations const &conv_sym,
                    double symprec) {
  std::optional<int> const aperiodic_axis = conv_prim.aperiodic_axis();
  ExactPositions out;
  out.reserve(static_cast<std::size_t>(conv_prim.size()));
  std::vector<int> indep_atoms;

  for (int i = 0; i < conv_prim.size(); ++i) {
    if (auto const eq = find_equivalent(i, indep_atoms, out, conv_prim,
                                        conv_sym, symprec)) {
      out.push_back({.position = eq->position,
                     .equivalent_atom = eq->representative});
    } else {
      indep_atoms.push_back(i);
      out.push_back({.position = set_exact_location(conv_prim.position(i),
                                                    conv_sym,
                                                    conv_prim.lattice(),
                                                    symprec, aperiodic_axis),
                     .equivalent_atom = i});
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
[[nodiscard]] std::optional<WyckoffLabel>
get_wyckoff_notation(Vector3d const &position, SymmetryOperations const &conv_sym,
                     int ref_multiplicity, Matrix3d const &lattice,
                     int hall_number, double symprec,
                     std::optional<int> aperiodic_axis) {
  std::vector<Vector3d> orbit;
  orbit.reserve(conv_sym.size());
  std::ranges::transform(conv_sym, std::back_inserter(orbit),
                         [&](auto const &op) { return op.apply(position); });

  data::WyckoffRange const range = data::wyckoff_indices(hall_number);
  for (int wi = 0; wi < range.count; ++wi) {
    data::WyckoffCoordinate const wc =
        data::wyckoff_coordinate(range.start + wi);

    for (auto const &oj : orbit) {
      int num_sitesym = 0;
      for (auto const &ok : orbit) {
        if (is_overlap(oj, ok, lattice, symprec, aperiodic_axis)) {
          Vector3d const mapped = wc.rotation.cast<double>() * ok + wc.translation;
          if (is_overlap(ok, mapped, lattice, symprec, aperiodic_axis)) {
            ++num_sitesym;
          }
        }
      }
      if (num_sitesym * wc.multiplicity == static_cast<int>(conv_sym.size()) &&
          wc.multiplicity == ref_multiplicity) {
        // The database stores Wyckoff positions in reverse order (g f e ... a).
        return WyckoffLabel{
            range.count - wi - 1,
            std::string(data::site_symmetry_symbol(range.start + wi))};
      }
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
