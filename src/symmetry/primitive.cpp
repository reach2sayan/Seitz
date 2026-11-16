#include "symmetry/primitive.hpp"

#include "symmetry/search.hpp"

#include "core/matrix_order.hpp"
#include "core/overlap.hpp"
#include "core/position_index.hpp"
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include <seitz/core/operation_set.hpp>
#include <seitz/core/periodicity.hpp>

#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

// The primitive cell, plus the cell-trimming helpers. Collect the pure
// translations; with only the trivial one the cell is already primitive and
// just needs Delaunay reduction. Otherwise pick three translation/unit vectors
// spanning the primitive volume, reduce, then fold and dedupe the atoms.
namespace seitz::symmetry {

namespace {

constexpr int kNumAttempt = 20;
constexpr double kReduceRate = 0.95;
constexpr int kTrimNumAttempt = 100;
constexpr double kTrimReduceRate = 0.95;
constexpr double kTrimIncreaseRate = 2.0;

[[nodiscard]] Vector3d row(Positions const &p, int i) {
  return p.row(i).transpose();
}

// Atom positions re-expressed through `to_new` and folded into [0, 1) on the
// periodic axes (the aperiodic axis of a layer cell is left raw).
[[nodiscard]] Positions transformed_positions(Matrix3d const &to_new,
                                              Cell const &cell) {
  Positions pos(cell.size(), 3);
  for (Index i = 0; i < cell.size(); ++i) {
    pos.row(i) = wrap(Vector3d(to_new * cell.position(i)), cell.periodicity())
                     .transpose();
  }
  return pos;
}

// The Delaunay reduction appropriate to the cell: 2D (periodic plane only) for
// a layer cell, full 3D otherwise.
template <GroupFamily F>
[[nodiscard]] Result<Lattice> reduce_lattice(Lattice const &lattice,
                                             CellPeriodicity const &periodicity,
                                             double symprec) {
  if constexpr (F == GroupFamily::layer) {
    return lattice.delaunay_in_plane(aperiodic_axis(periodicity).value_or(2),
                                     symprec);
  } else {
    return lattice.delaunay(symprec);
  }
}

// The multiplicity-one case: reduce the lattice and fold the atoms into it.
template <GroupFamily F>
[[nodiscard]] std::optional<Cell> smallest_lattice_cell(Cell const &cell,
                                                        double symprec) {
  auto const min_lat =
      reduce_lattice<F>(cell.lattice(), cell.periodicity(), symprec);
  if (!min_lat) {
    return std::nullopt;
  }
  Matrix3d const trans = min_lat->matrix().inverse() * cell.lattice().matrix();
  return Cell(*min_lat, transformed_positions(trans, cell), cell.types(),
              cell.periodicity());
}

// Clean a candidate relative lattice (whose |det| should equal `multi`) via its
// exact integer inverse, then Delaunay-reduce and return the primitive lattice.
// Shared by the 3D and layer paths.
template <GroupFamily F>
[[nodiscard]] std::optional<Lattice>
finish_primitive_lattice(Matrix3d relative, Lattice const &cell_lattice,
                         int multi, CellPeriodicity const &periodicity,
                         double symprec) {
  auto const rel_inv = math::inverse(relative, symprec);
  if (rel_inv) {
    Matrix3i const inv_int = math::round_to_int(*rel_inv);
    if (std::abs(inv_int.determinant()) == multi) {
      relative = inv_int.cast<double>().inverse();
    }
  }
  auto const reduced = reduce_lattice<F>(cell_lattice.transformed(relative),
                                         periodicity, symprec);
  return reduced.has_value() ? std::optional(reduced.value())
                             : std::nullopt;
}

template <GroupFamily F>
[[nodiscard]] std::optional<Lattice>
primitive_lattice(Cell const &cell, std::span<Vector3d const> pure_trans,
                  double symprec) {
  int const multi = static_cast<int>(pure_trans.size());
  double const init_volume = cell.lattice().volume();
  // Layer: the aperiodic lattice vector keeps its column so the 2D reduction
  // leaves it alone; the periodic ones come from the candidates below.
  std::optional<int> const ap = [&]() -> std::optional<int> {
    if constexpr (F == GroupFamily::layer) {
      return aperiodic_axis(cell.periodicity()).value_or(2);
    } else {
      return std::nullopt;
    }
  }();
  // Candidate basis vectors: the pure translations (the zero translation, the
  // reference atom's own, first) and the periodic unit vectors.
  std::vector<Vector3d> cand(pure_trans.begin(), pure_trans.end());
  for (int a = 0; a < 3; ++a) {
    if (a != ap) {
      cand.push_back(Vector3d::Unit(a));
    }
  }
  auto const right_volume = [&](Matrix3d const &relative) {
    double const volume =
        std::abs((cell.lattice().matrix() * relative).determinant());
    return volume > symprec && math::nint(init_volume / volume) == multi;
  };
  auto const finish = [&](Matrix3d const &relative) {
    return finish_primitive_lattice<F>(relative, cell.lattice(), multi,
                                       cell.periodicity(), symprec);
  };

  if constexpr (F == GroupFamily::layer) {
    std::array<int, 2> const periodic{*ap == 0 ? 1 : 0, *ap == 2 ? 1 : 2};
    auto const ids = std::views::iota(std::size_t{0}, cand.size());
    // The first pair (in lexicographic order) spanning a cell of the right
    // volume.
    for (auto const [i, j] : std::views::cartesian_product(ids, ids) |
                                 std::views::filter([](auto const &ij) {
                                   return std::get<0>(ij) < std::get<1>(ij);
                                 })) {
      Matrix3d relative;
      relative.col(*ap) = Vector3d::Unit(*ap);
      relative.col(periodic[0]) = cand[i];
      relative.col(periodic[1]) = cand[j];
      if (right_volume(relative)) {
        return finish(relative);
      }
    }
  } else {
    // The first triple (in lexicographic order) of non-zero candidates
    // spanning a cell of the right volume. The zero candidate gives an exactly
    // zero determinant, so it is filtered out before the product rather than
    // visited in the O(n^2) triples that hold it.
    std::vector<std::size_t> const nonzero{
        std::from_range, std::views::iota(std::size_t{0}, cand.size()) |
                             std::views::filter([&](std::size_t i) {
                               return !(cand[i].array() == 0.0).all();
                             })};
    for (auto const [i, j, k] :
         std::views::cartesian_product(nonzero, nonzero, nonzero) |
             std::views::filter([](auto const &ijk) {
               auto const [i, j, k] = ijk;
               return i < j && j < k;
             })) {
      Matrix3d relative;
      relative << cand[i], cand[j], cand[k];
      if (right_volume(relative)) {
        return finish(relative);
      }
    }
  }
  return std::nullopt;
}

// Fold atoms into the trimmed lattice, de-dup translationally-equivalent ones,
// average positions.
[[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
trim_cell(Lattice const &trimmed_lattice, Cell const &cell, double symprec) {
  int const n = static_cast<int>(cell.size());
  int const ratio =
      std::abs(math::nint(cell.lattice().matrix().determinant() /
                          trimmed_lattice.matrix().determinant()));
  if (ratio == 0 || n % ratio != 0) {
    return std::nullopt;
  }

  Matrix3i const tmat = math::round_to_int(
      Matrix3d(trimmed_lattice.matrix().inverse() * cell.lattice().matrix()));
  if (std::abs(tmat.determinant()) != ratio) {
    return std::nullopt;
  }

  // Atom positions expressed in the trimmed basis, folded into [0, 1) on the
  // periodic axes (the aperiodic axis is left raw for a layer cell).
  Positions const pos = transformed_positions(tmat.cast<double>(), cell);
  CellPeriodicity const &periodicity = cell.periodicity();

  // Representatives, with tolerance adjusted until each class has `ratio`
  // atoms: overlap[i] is the lowest-index representative of i's type coinciding
  // with i, or i itself, which then starts a class.
  std::vector<int> overlap(static_cast<std::size_t>(n));
  // Hoisted out of the retry loop below: the keys change every attempt but the
  // capacity does not, and this loop runs up to kTrimNumAttempt times.
  boost::container::flat_map<int, int> class_size;
  class_size.reserve(static_cast<std::size_t>(n));
  double tol = symprec;
  bool ok = false;
  PositionIndex::Scratch scratch;
  for (int attempt = 0; attempt < kTrimNumAttempt && !ok; ++attempt) {
    PositionIndex const index(pos, cell.types(), trimmed_lattice.matrix(), tol,
                              periodicity);
    // Representatives claim, in index order, every unclaimed later atom of
    // their type they coincide with; an atom still unclaimed when its turn
    // comes is the next representative. Same classes as asking each atom for
    // its lowest coinciding representative (coincidence is symmetric), but
    // one query per representative rather than per atom -- a supercell folds
    // `ratio` atoms onto each site, and each of those queries returned them
    // all.
    std::ranges::iota(overlap, 0);
    auto const unclaimed = [&](int i) {
      return overlap[static_cast<std::size_t>(i)] == i;
    };
    for (int const i : std::views::iota(0, n) | std::views::filter(unclaimed)) {
      for (int const k : index.matches(row(pos, i), cell.type(i), scratch) |
                             std::views::filter([&](int k) {
                               return k > i && unclaimed(k);
                             })) {
        overlap[static_cast<std::size_t>(k)] = i;
      }
    }

    // Class sizes in one pass: every overlap entry names its representative,
    // so the histogram keys are exactly the representatives, in index order.
    class_size.clear();
    for (int const rep : overlap) {
      ++class_size[rep];
    }
    auto const bad = std::ranges::find_if(
        class_size, [&](auto const &entry) { return entry.second != ratio; });
    ok = bad == class_size.end();
    if (!ok) {
      // The smallest-index wrong-sized class decides the adjustment: too few
      // atoms widens the tolerance, too many tightens it.
      tol *= bad->second < ratio ? kTrimIncreaseRate : kTrimReduceRate;
    }
  }
  if (!ok) {
    return std::nullopt;
  }

  // Build the mapping (input atom -> trimmed atom) and the trimmed types:
  // each representative becomes the next trimmed atom, every other atom maps
  // where its representative went.
  std::vector<int> mapping(static_cast<std::size_t>(n));
  std::vector<int> trimmed_types;
  for (auto const [i, rep] : overlap | std::views::enumerate) {
    auto const ui = static_cast<std::size_t>(i);
    if (rep == i) {
      mapping[ui] = static_cast<int>(trimmed_types.size());
      trimmed_types.push_back(cell.type(i));
    } else {
      mapping[ui] = mapping[static_cast<std::size_t>(rep)];
    }
  }
  int const tn = static_cast<int>(trimmed_types.size());

  // Average the positions of overlapping atoms (with periodic-boundary care):
  // a component more than half a cell from its representative is shifted one
  // cell toward it before summing.
  Positions tpos = Positions::Zero(tn, 3);
  for (auto const [i, j, k] :
       std::views::zip(std::views::iota(0, n), mapping, overlap)) {
    tpos.row(j) += pos.row(i).binaryExpr(pos.row(k), [](double pi, double pk) {
      if (std::abs(pk - pi) > 0.5) {
        return pi < pk ? pi + 1.0 : pi - 1.0;
      }
      return pi;
    });
  }
  int const multi = n / tn;
  for (int i = 0; i < tn; ++i) {
    tpos.row(i) = wrap(Vector3d(tpos.row(i).transpose() / multi), periodicity)
                      .transpose();
  }

  return std::pair{Cell(trimmed_lattice, std::move(tpos),
                        std::move(trimmed_types), periodicity),
                   std::move(mapping)};
}

// The primitive lattice in "translation space": reduce a cell whose atoms sit
// at the pure translations, and its lattice is the
// (primitive-to-conventional)^-1 transformation. nullopt unless they span
// exactly one primitive point.
[[nodiscard]] std::optional<Matrix3d>
primitive_in_translation_space(std::vector<Vector3d> const &pure_trans,
                               std::size_t symmetry_size, double symprec) {
  std::size_t const np = pure_trans.size();
  if (np == 0 || symmetry_size % np != 0) {
    return std::nullopt;
  }
  Cell const cell(Lattice{Matrix3d::Identity()}, to_positions(pure_trans),
                  Types(np, 1));
  auto const prim =
      PrimitiveFinder<GroupFamily::space>{cell, {symprec, std::nullopt}}.find();
  if (!prim || prim->cell.size() != 1) {
    return std::nullopt;
  }
  return prim->cell.lattice().matrix();
}

// The first occurrence of each distinct rotation, keeping its translation.
// std::nullopt if the number of distinct rotations differs from `primsym_size`.
[[nodiscard]] std::optional<Operations>
collect_primitive_operations(std::span<SymmetryOperation const> operations,
                             std::size_t primsym_size) {
  auto out = unique_by_rotation(operations, &SymmetryOperation::rotation);
  if (out.size() != primsym_size) {
    return std::nullopt;
  }
  return Operations{std::move(out)};
}

} // namespace

template <GroupFamily F>
std::optional<std::pair<Cell, std::vector<int>>>
PrimitiveFinder<F>::trim_to(Lattice const &trimmed_lattice) const {
  return trim_cell(trimmed_lattice, cell_, tol_.symprec);
}

// The primitive cell implied by `pure_trans` at one tolerance: a single pure
// translation means the cell is already primitive and only needs reducing.
template <GroupFamily F>
[[nodiscard]] std::optional<Primitive>
primitive_from(Cell const &cell, std::span<Vector3d const> pure_trans,
               Tolerance const &tol) {
  if (pure_trans.size() == 1) {
    auto smallest = smallest_lattice_cell<F>(cell, tol.symprec);
    if (!smallest) {
      return std::nullopt;
    }
    std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
    std::ranges::iota(mapping, 0);
    return Primitive{std::move(*smallest), std::move(mapping),
                     cell.lattice().matrix(), tol};
  }
  auto const prim_lat = primitive_lattice<F>(cell, pure_trans, tol.symprec);
  if (!prim_lat) {
    return std::nullopt;
  }
  auto trimmed = trim_cell(*prim_lat, cell, tol.symprec);
  if (!trimmed) {
    return std::nullopt;
  }
  return Primitive{std::move(trimmed->first), std::move(trimmed->second),
                   cell.lattice().matrix(), tol};
}

template <GroupFamily F> Result<Primitive> PrimitiveFinder<F>::find() const {
  double tolerance = tol_.symprec;
  for (int attempt = 0; attempt < kNumAttempt;
       ++attempt, tolerance *= kReduceRate) {
    Tolerance const tol{tolerance, tol_.angle_tolerance};
    auto const pure = SymmetrySearch<F>{cell_, tol}.pure_translations();
    if (pure.empty()) {
      continue;
    }
    if (auto prim = primitive_from<F>(cell_, pure, tol)) {
      return std::move(*prim);
    }
  }
  return leaf::new_error(e_cell_standardization_failed{});
}

template <GroupFamily F>
std::optional<Lattice> PrimitiveFinder<F>::lattice_from_pure_translations(
    std::span<Vector3d const> pure_trans) const {
  return primitive_lattice<F>(cell_, pure_trans, tol_.symprec);
}

template <GroupFamily F>
Result<Primitive> PrimitiveFinder<F>::from_pure_translations(
    std::span<Vector3d const> pure_trans) const {
  if (auto prim = primitive_from<F>(cell_, pure_trans,
                                    {tol_.symprec, std::nullopt})) {
    return std::move(*prim);
  }
  return leaf::new_error(e_cell_standardization_failed{});
}

template class PrimitiveFinder<GroupFamily::space>;
template class PrimitiveFinder<GroupFamily::layer>;

} // namespace seitz::symmetry

// The implementation behind OperationSet::to_primitive, out of line here so
// the public header never names the primitive-cell machinery.
namespace seitz::detail {

std::optional<std::pair<Operations, Matrix3d>>
primitive_operations(Operations const &operations, Tolerance const &tol) {
  auto const pure_trans = operations.pure_translations();
  if (pure_trans.empty()) {
    return std::nullopt;
  }
  std::size_t const primsym_size = operations.size() / pure_trans.size();

  // t_mat transforms primitive -> conventional; t_mat_inv is its inverse,
  // recovered as the primitive lattice in translation space.
  auto const t_mat_inv = symmetry::primitive_in_translation_space(
      pure_trans, operations.size(), tol.symprec);
  if (!t_mat_inv) {
    return std::nullopt;
  }
  Matrix3d const t_mat = t_mat_inv->inverse();

  auto const prim =
      symmetry::collect_primitive_operations(operations.span(), primsym_size);
  if (!prim) {
    return std::nullopt;
  }

  // (T, 0) (R, t) (T, 0)^-1 = (T R T^-1, T t).
  return std::pair{prim->conjugated_by(t_mat, *t_mat_inv), t_mat};
}

} // namespace seitz::detail
