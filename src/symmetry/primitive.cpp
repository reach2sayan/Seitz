#include <cppcrystal/symmetry/primitive.hpp>

#include <cppcrystal/core/matrix_order.hpp>
#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/core/position_index.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

// Find the primitive cell, plus the cell-trimming helpers:
//   1. Collect the pure translations of the cell.
//   2. If there is only the trivial one, the cell is already primitive; just
//      Delaunay-reduce its lattice.
//   3. Otherwise pick three pure-translation/unit vectors that span the
//      primitive volume, Delaunay-reduce, then fold and de-duplicate the atoms
//      into the smaller cell.
namespace cppcrystal::symmetry {

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
[[nodiscard]] Result<Lattice> reduce_lattice(Lattice const &lattice,
                                             CellPeriodicity const &periodicity,
                                             double symprec) {
  auto const axis = aperiodic_axis(periodicity);
  return axis ? lattice.delaunay(*axis, symprec) : lattice.delaunay(symprec);
}

// The multiplicity-one case: reduce the lattice and fold the atoms into it.
[[nodiscard]] std::optional<Cell> smallest_lattice_cell(Cell const &cell,
                                                        double symprec) {
  auto const min_lat = reduce_lattice(cell.lattice(), cell.periodicity(), symprec);
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
  auto const reduced =
      reduce_lattice(cell_lattice.transformed(relative), periodicity, symprec);
  return reduced.has_value() ? std::optional<Lattice>(reduced.value())
                             : std::nullopt;
}

[[nodiscard]] std::optional<Lattice>
primitive_lattice(Cell const &cell, std::vector<Vector3d> const &pure_trans,
                  double symprec) {
  int const multi = static_cast<int>(pure_trans.size());
  double const init_volume = cell.lattice().volume();

  if (auto const layer_axis = aperiodic_axis(cell.periodicity())) {
    // Layer: the third basis vector is fixed to the aperiodic lattice vector
    // (the cell is not periodic along it); the two periodic vectors are chosen
    // from the in-plane pure translations and the other two unit vectors. The
    // aperiodic axis is kept at its own column index so the 2D reduction below
    // leaves it untouched.
    int const ap = *layer_axis;
    std::vector<Vector3d> cand = pure_trans; // in-plane, includes zero
    for (int a = 0; a < 3; ++a) {
      if (a != ap) {
        cand.push_back(Vector3d::Unit(a));
      }
    }
    std::array<int, 2> const periodic{ap == 0 ? 1 : 0, ap == 2 ? 1 : 2};
    auto const ids = std::views::iota(std::size_t{0}, cand.size());
    // The first pair (in lexicographic order) spanning a cell of the right
    // volume.
    for (auto const [i, j] : std::views::cartesian_product(ids, ids) |
                                 std::views::filter([](auto const &ij) {
                                   return std::get<0>(ij) < std::get<1>(ij);
                                 })) {
      Matrix3d relative;
      relative.col(ap) = Vector3d::Unit(ap);
      relative.col(periodic[0]) = cand[i];
      relative.col(periodic[1]) = cand[j];
      double const volume = std::abs((cell.lattice().matrix() * relative).determinant());
      if (volume <= symprec || math::nint(init_volume / volume) != multi) {
        continue;
      }
      return finish_primitive_lattice(relative, cell.lattice(), multi,
                                      cell.periodicity(), symprec);
    }
    return std::nullopt;
  }

  std::vector<Vector3d> cand = pure_trans; // includes the zero translation
  cand.push_back(Vector3d::UnitX());
  cand.push_back(Vector3d::UnitY());
  cand.push_back(Vector3d::UnitZ());
  auto const ids = std::views::iota(std::size_t{0}, cand.size());
  // The first triple (in lexicographic order) spanning a cell of the right
  // volume.
  for (auto const [i, j, k] :
       std::views::cartesian_product(ids, ids, ids) |
           std::views::filter([](auto const &ijk) {
             return std::get<0>(ijk) < std::get<1>(ijk) &&
                    std::get<1>(ijk) < std::get<2>(ijk);
           })) {
    Matrix3d relative;
    relative.col(0) = cand[i];
    relative.col(1) = cand[j];
    relative.col(2) = cand[k];
    double const volume = std::abs((cell.lattice().matrix() * relative).determinant());
    if (volume <= symprec || math::nint(init_volume / volume) != multi) {
      continue;
    }
    return finish_primitive_lattice(relative, cell.lattice(), multi,
                                    cell.periodicity(), symprec);
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

  // Representative table with tolerance adjustment until each class has
  // `ratio` atoms: overlap[i] is the lowest-index class representative (an
  // atom that is its own representative) of i's type coinciding with i, or i
  // itself, which then starts a class.
  std::vector<int> overlap(static_cast<std::size_t>(n));
  double tol = symprec;
  bool ok = false;
  for (int attempt = 0; attempt < kTrimNumAttempt && !ok; ++attempt) {
    PositionIndex const index(
        BucketGeometry::of(trimmed_lattice.matrix(), tol, periodicity), pos,
        cell.types(), trimmed_lattice.matrix(), tol, periodicity);
    for (int i = 0; i < n; ++i) {
      auto const ui = static_cast<std::size_t>(i);
      overlap[ui] = i; // i is a representative until a lower one claims it
      overlap[ui] = index
                        .first_match(row(pos, i), cell.type(i),
                                     [&](int j) {
                                       return j <= i &&
                                              overlap[static_cast<std::size_t>(
                                                  j)] == j;
                                     })
                        .value_or(i);
    }

    // Class sizes in one pass: every overlap entry names its representative,
    // so the histogram keys are exactly the representatives, in index order.
    std::map<int, int> class_size;
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
  if (!ok)
    return std::nullopt;

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
  for (int i = 0; i < n; ++i) {
    int const j = mapping[static_cast<std::size_t>(i)];
    int const k = overlap[static_cast<std::size_t>(i)];
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

  return std::make_pair(
      Cell(trimmed_lattice, tpos, trimmed_types, periodicity), mapping);
}

// The translations of a symmetry-operation set whose rotation is the identity.
[[nodiscard]] std::vector<Vector3d>
operation_pure_translations(SymmetryOperations const &operations) {
  std::vector<Vector3d> out;
  std::ranges::copy(operations | std::views::filter([](auto const &op) {
                      return op.rotation == Matrix3i::Identity();
                    }) | std::views::transform(&SymmetryOperation::translation),
                    std::back_inserter(out));
  return out;
}

// The primitive lattice in "translation space": a unit cell whose atoms sit at
// the pure translations is reduced to its primitive cell; that cell's lattice
// is the (primitive-to-conventional)^-1 transformation. std::nullopt unless the
// translations span exactly one primitive point.
[[nodiscard]] std::optional<Matrix3d>
primitive_in_translation_space(std::vector<Vector3d> const &pure_trans,
                               std::size_t symmetry_size, double symprec) {
  std::size_t const np = pure_trans.size();
  if (np == 0 || symmetry_size % np != 0) {
    return std::nullopt;
  }
  Cell const cell(Lattice{Matrix3d::Identity()}, to_positions(pure_trans),
                  Types(np, 1));
  auto const prim = find_primitive(cell, Tolerance{symprec, std::nullopt});
  if (!prim || prim->cell.size() != 1) {
    return std::nullopt;
  }
  return prim->cell.lattice().matrix();
}

// The first occurrence of each distinct rotation, keeping its translation.
// std::nullopt if the number of distinct rotations differs from `primsym_size`.
[[nodiscard]] std::optional<SymmetryOperations>
collect_primitive_operations(SymmetryOperations const &operations,
                             std::size_t primsym_size) {
  auto out = unique_by_rotation(operations, &SymmetryOperation::rotation);
  if (out.size() != primsym_size) {
    return std::nullopt;
  }
  return out;
}

} // namespace

std::optional<std::pair<Cell, std::vector<int>>>
trim_to_lattice(Lattice const &trimmed_lattice, Cell const &cell,
                double symprec) {
  return trim_cell(trimmed_lattice, cell, symprec);
}

Result<Primitive> find_primitive(Cell const &cell, Tolerance const &tol) {
  double tolerance = tol.symprec;
  for (int attempt = 0; attempt < kNumAttempt;
       ++attempt, tolerance *= kReduceRate) {
    auto const pure = pure_translations(cell, tolerance);
    if (pure.empty()) {
      continue;
    }

    if (pure.size() == 1) {
      auto smallest = smallest_lattice_cell(cell, tolerance);
      if (!smallest) {
        continue;
      }
      std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
      std::ranges::iota(mapping, 0);
      return Primitive{std::move(*smallest), std::move(mapping),
                       cell.lattice().matrix(),
                       {tolerance, tol.angle_tolerance}};
    }

    auto const prim_lat = primitive_lattice(cell, pure, tolerance);
    if (!prim_lat) {
      continue;
    }
    auto trimmed = trim_cell(*prim_lat, cell, tolerance);
    if (!trimmed) {
      continue;
    }
    return Primitive{std::move(trimmed->first), std::move(trimmed->second),
                     cell.lattice().matrix(),
                     {tolerance, tol.angle_tolerance}};
  }
  return leaf::new_error(e_cell_standardization_failed{});
}

std::optional<Lattice> primitive_lattice_vectors(
    Cell const &cell, std::vector<Vector3d> const &pure_trans, double symprec) {
  return primitive_lattice(cell, pure_trans, symprec);
}

std::optional<std::pair<SymmetryOperations, Matrix3d>>
primitive_symmetry(SymmetryOperations const &operations, double symprec) {
  auto const pure_trans = operation_pure_translations(operations);
  if (pure_trans.empty()) {
    return std::nullopt;
  }
  std::size_t const primsym_size = operations.size() / pure_trans.size();

  // t_mat transforms primitive -> conventional; t_mat_inv is its inverse,
  // recovered as the primitive lattice in translation space.
  auto const t_mat_inv =
      primitive_in_translation_space(pure_trans, operations.size(), symprec);
  if (!t_mat_inv) {
    return std::nullopt;
  }
  Matrix3d const t_mat = t_mat_inv->inverse();

  auto prim = collect_primitive_operations(operations, primsym_size);
  if (!prim) {
    return std::nullopt;
  }

  // (T, 0) (R, t) (T, 0)^-1 = (T R T^-1, T t).
  for (auto &op : *prim) {
    op = conjugated_by(op, t_mat, *t_mat_inv);
  }
  return std::make_pair(std::move(*prim), t_mat);
}

Result<Primitive>
find_primitive_with_pure_translations(Cell const &cell,
                                      std::vector<Vector3d> const &pure_trans,
                                      double symprec) {
  if (pure_trans.size() == 1) {
    auto smallest = smallest_lattice_cell(cell, symprec);
    if (!smallest) {
      return leaf::new_error(e_cell_standardization_failed{});
    }
    std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
    std::iota(mapping.begin(), mapping.end(), 0);
    return Primitive{std::move(*smallest), std::move(mapping),
                     cell.lattice().matrix(), {symprec, std::nullopt}};
  }
  auto const prim_lat = primitive_lattice(cell, pure_trans, symprec);
  if (!prim_lat) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  auto trimmed = trim_cell(*prim_lat, cell, symprec);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return Primitive{std::move(trimmed->first), std::move(trimmed->second),
                   cell.lattice().matrix(), {symprec, std::nullopt}};
}

} // namespace cppcrystal::symmetry
