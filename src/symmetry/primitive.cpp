#include <spglib/symmetry/primitive.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>
#include <spglib/reduce/delaunay.hpp>
#include <spglib/symmetry/find_symmetry.hpp>

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

// Port of primitive.c (3D space-group path) + the cell-trimming helpers from
// cell.c (trim_cell / get_overlap_table / translate_atoms_in_trimmed_lattice).
//   1. Collect the pure translations of the cell.
//   2. If there is only the trivial one, the cell is already primitive; just
//      Delaunay-reduce its lattice.
//   3. Otherwise pick three pure-translation/unit vectors that span the
//      primitive volume, Delaunay-reduce, then fold and de-duplicate the atoms
//      into the smaller cell.
namespace spglib::symmetry {

namespace {

constexpr int kNumAttempt = 20; // primitive.c get_primitive
constexpr double kReduceRate = 0.95;
constexpr int kTrimNumAttempt = 100; // cell.c trim
constexpr double kTrimReduceRate = 0.95;
constexpr double kTrimIncreaseRate = 2.0;

[[nodiscard]] Vector3d row(Positions const &p, int i) {
  return p.row(i).transpose();
}

// Fold a fractional coordinate into the cell, but leave the aperiodic axis (if
// any) at its raw value — layer cells are not periodic along it.
[[nodiscard]] Vector3d wrap_periodic(Vector3d const &v,
                                     std::optional<int> aperiodic_axis) {
  Vector3d out;
  for (int j = 0; j < 3; ++j) {
    out[j] = (aperiodic_axis && j == *aperiodic_axis)
                 ? v[j]
                 : math::wrap_to_unit_cell(v[j]);
  }
  return out;
}

// The Delaunay reduction appropriate to the cell: 2D (periodic plane only) for
// a layer cell, full 3D otherwise.
[[nodiscard]] Result<Matrix3d>
reduce_lattice(Matrix3d const &lattice, std::optional<int> aperiodic_axis,
               double symprec) {
  return aperiodic_axis
             ? reduce::delaunay_reduce(lattice, *aperiodic_axis, symprec)
             : reduce::delaunay_reduce(lattice, symprec);
}

// get_cell_with_smallest_lattice: the multiplicity-one case.
[[nodiscard]] std::optional<Cell>
smallest_lattice_cell(Cell const &cell, std::optional<int> aperiodic_axis,
                      double symprec) {
  auto const min_lat = reduce_lattice(cell.lattice(), aperiodic_axis, symprec);
  if (!min_lat) {
    return std::nullopt;
  }
  Matrix3d const trans = min_lat->inverse() * cell.lattice();
  Positions pos(cell.size(), 3);
  for (Index i = 0; i < cell.size(); ++i) {
    pos.row(i) =
        wrap_periodic(Vector3d(trans * cell.position(i)), aperiodic_axis)
            .transpose();
  }
  return Cell(*min_lat, pos, cell.types(), aperiodic_axis);
}

// find_primitive_lattice_vectors + cleaning + Delaunay reduce. Picks the first
// triple of {pure translations, unit vectors} spanning the primitive volume.
// Clean a candidate relative lattice (whose |det| should equal `multi`) via its
// exact integer inverse, then reduce and return the primitive lattice. Shared by
// the 3D and layer paths.
[[nodiscard]] std::optional<Matrix3d>
finish_primitive_lattice(Matrix3d relative, Matrix3d const &cell_lattice,
                         int multi, std::optional<int> aperiodic_axis,
                         double symprec) {
  auto const rel_inv = math::inverse(relative, symprec);
  if (rel_inv) {
    Matrix3i const inv_int = math::round_to_int(*rel_inv);
    if (std::abs(inv_int.determinant()) == multi) {
      relative = inv_int.cast<double>().inverse();
    }
  }
  Matrix3d const prim = cell_lattice * relative;
  auto const reduced = reduce_lattice(prim, aperiodic_axis, symprec);
  return reduced.has_value() ? std::optional<Matrix3d>(reduced.value())
                             : std::nullopt;
}

[[nodiscard]] std::optional<Matrix3d>
primitive_lattice(Cell const &cell, std::vector<Vector3d> const &pure_trans,
                  std::optional<int> aperiodic_axis, double symprec) {
  int const multi = static_cast<int>(pure_trans.size());
  double const init_volume = std::abs(cell.lattice().determinant());

  if (aperiodic_axis) {
    // Layer: the third basis vector is fixed to the aperiodic lattice vector
    // (the cell is not periodic along it); the two periodic vectors are chosen
    // from the in-plane pure translations and the other two unit vectors. The
    // aperiodic axis is kept at its own column index so the 2D reduction below
    // leaves it untouched. (primitive.c find_primitive_lattice_vectors, layer.)
    int const ap = *aperiodic_axis;
    std::vector<Vector3d> cand = pure_trans; // in-plane, includes zero
    for (int a = 0; a < 3; ++a) {
      if (a != ap) {
        cand.push_back(Vector3d::Unit(a));
      }
    }
    std::array<int, 2> const periodic{ap == 0 ? 1 : 0, ap == 2 ? 1 : 2};
    std::size_t const n = cand.size();
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        Matrix3d relative;
        relative.col(ap) = Vector3d::Unit(ap);
        relative.col(periodic[0]) = cand[i];
        relative.col(periodic[1]) = cand[j];
        double const volume =
            std::abs((cell.lattice() * relative).determinant());
        if (volume <= symprec || math::nint(init_volume / volume) != multi) {
          continue;
        }
        return finish_primitive_lattice(relative, cell.lattice(), multi,
                                        aperiodic_axis, symprec);
      }
    }
    return std::nullopt;
  }

  std::vector<Vector3d> cand = pure_trans; // includes the zero translation
  cand.push_back(Vector3d::UnitX());
  cand.push_back(Vector3d::UnitY());
  cand.push_back(Vector3d::UnitZ());
  std::size_t const n = cand.size();

  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j)
      for (std::size_t k = j + 1; k < n; ++k) {
        Matrix3d tmp;
        tmp.col(0) = cell.lattice() * cand[i];
        tmp.col(1) = cell.lattice() * cand[j];
        tmp.col(2) = cell.lattice() * cand[k];
        double const volume = std::abs(tmp.determinant());
        if (volume <= symprec || math::nint(init_volume / volume) != multi) {
          continue;
        }

        Matrix3d relative;
        relative.col(0) = cand[i];
        relative.col(1) = cand[j];
        relative.col(2) = cand[k];
        // matches spglib: first valid triple only
        return finish_primitive_lattice(relative, cell.lattice(), multi,
                                        std::nullopt, symprec);
      }
  return std::nullopt;
}

// trim_cell + get_overlap_table: fold atoms into the trimmed lattice, de-dup
// translationally-equivalent ones, average positions.
[[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
trim_cell(Matrix3d const &trimmed_lattice, Cell const &cell,
          std::optional<int> aperiodic_axis, double symprec) {
  int const n = static_cast<int>(cell.size());
  int const ratio = std::abs(
      math::nint(cell.lattice().determinant() / trimmed_lattice.determinant()));
  if (ratio == 0 || n % ratio != 0) {
    return std::nullopt;
  }

  Matrix3i const tmat =
      math::round_to_int(Matrix3d(trimmed_lattice.inverse() * cell.lattice()));
  if (std::abs(tmat.determinant()) != ratio) {
    return std::nullopt;
  }

  // Atom positions expressed in the trimmed basis, folded into [0, 1) on the
  // periodic axes (the aperiodic axis is left raw for a layer cell).
  Positions pos(n, 3);
  for (int i = 0; i < n; ++i) {
    pos.row(i) = wrap_periodic(Vector3d(tmat.cast<double>() * cell.position(i)),
                               aperiodic_axis)
                     .transpose();
  }

  // Overlap table with tolerance adjustment until each class has `ratio` atoms.
  std::vector<int> overlap(static_cast<std::size_t>(n));
  double tol = symprec;
  bool ok = false;
  for (int attempt = 0; attempt < kTrimNumAttempt && !ok; ++attempt) {

    for (int i = 0; i < n; ++i) {
      overlap[static_cast<std::size_t>(i)] = i;
      for (int j = 0; j < n; ++j) {
        if (cell.type(i) == cell.type(j) &&
            overlap[static_cast<std::size_t>(j)] == j &&
            is_overlap(row(pos, i), row(pos, j), trimmed_lattice, tol,
                       aperiodic_axis)) {
          overlap[static_cast<std::size_t>(i)] = j;
          break;
        }
      }
    }

    bool retry = false;
    for (int i = 0; i < n && !retry; ++i) {
      if (overlap[static_cast<std::size_t>(i)] != i) {
        continue;
      }
      auto const count = std::ranges::count_if(
          std::views::iota(0, n),
          [&](int j) { return overlap[static_cast<std::size_t>(j)] == i; });
      if (count < ratio) {
        tol *= kTrimIncreaseRate;
        retry = true;
      } else if (count > ratio) {
        tol *= kTrimReduceRate;
        retry = true;
      }
    }
    ok = !retry;
  }
  if (!ok)
    return std::nullopt;

  // Build the mapping (input atom -> trimmed atom) and the trimmed types.
  std::vector<int> mapping(static_cast<std::size_t>(n));
  std::vector<int> trimmed_types;
  int index_atom = 0;
  for (int i = 0; i < n; ++i) {
    auto const ui = static_cast<std::size_t>(i);
    if (overlap[ui] == i) {
      mapping[ui] = index_atom++;
      trimmed_types.push_back(cell.type(i));
    } else {
      mapping[ui] = mapping[static_cast<std::size_t>(overlap[ui])];
    }
  }
  int const tn = index_atom;

  // Average the positions of overlapping atoms (with periodic-boundary care).
  Positions tpos = Positions::Zero(tn, 3);
  for (int i = 0; i < n; ++i) {
    int const j = mapping[static_cast<std::size_t>(i)];
    int const k = overlap[static_cast<std::size_t>(i)];
    for (int l = 0; l < 3; ++l) {
      double const pi = pos(i, l);
      double const pk = pos(k, l);
      if (std::abs(pk - pi) > 0.5) {
        tpos(j, l) += pi < pk ? pi + 1.0 : pi - 1.0;
      } else {
        tpos(j, l) += pi;
      }
    }
  }
  int const multi = n / tn;
  for (int i = 0; i < tn; ++i) {
    tpos.row(i) =
        wrap_periodic(Vector3d(tpos.row(i).transpose() / multi), aperiodic_axis)
            .transpose();
  }

  return std::make_pair(
      Cell(trimmed_lattice, tpos, trimmed_types, aperiodic_axis), mapping);
}

// The translations of a symmetry-operation set whose rotation is the identity
// (primitive.c collect_pure_translations, the symmetry-set variant).
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
// translations span exactly one primitive point. Port of
// get_primitive_in_translation_space.
[[nodiscard]] std::optional<Matrix3d>
primitive_in_translation_space(std::vector<Vector3d> const &pure_trans,
                               std::size_t symmetry_size, double symprec) {
  std::size_t const np = pure_trans.size();
  if (np == 0 || symmetry_size % np != 0) {
    return std::nullopt;
  }
  Positions pos(static_cast<Index>(np), 3);
  for (auto const [i, t] : pure_trans | std::views::enumerate) {
    pos.row(static_cast<Index>(i)) = t.transpose();
  }
  Cell const cell(Matrix3d::Identity(), pos, Types(np, 1));
  auto const prim = find_primitive(cell, symprec);
  if (!prim || prim->cell.size() != 1) {
    return std::nullopt;
  }
  return prim->cell.lattice();
}

// The first occurrence of each distinct rotation, keeping its translation,
// until `primsym_size` are collected. std::nullopt if the number of distinct
// rotations differs from `primsym_size`. Port of collect_primitive_symmetry.
[[nodiscard]] std::optional<SymmetryOperations>
collect_primitive_operations(SymmetryOperations const &operations,
                             std::size_t primsym_size) {
  SymmetryOperations out;
  for (auto const &op : operations) {
    if (std::ranges::any_of(out, [&](SymmetryOperation const &o) {
          return o.rotation == op.rotation;
        })) {
      continue;
    }
    if (out.size() == primsym_size) {
      return std::nullopt; // more distinct rotations than expected
    }
    out.push_back(op);
  }
  if (out.size() != primsym_size) {
    return std::nullopt;
  }
  return out;
}

} // namespace

Result<Primitive> find_primitive(Cell const &cell, double symprec,
                                 AngleTolerance angle_tolerance) {
  std::optional<int> const aperiodic_axis = cell.aperiodic_axis();
  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttempt;
       ++attempt, tolerance *= kReduceRate) {
    auto const pure = pure_translations(cell, tolerance);
    if (pure.empty()) {
      continue;
    }

    if (pure.size() == 1) {
      auto smallest = smallest_lattice_cell(cell, aperiodic_axis, tolerance);
      if (!smallest) {
        continue;
      }
      std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
      std::iota(mapping.begin(), mapping.end(), 0);
      return Primitive{std::move(*smallest), std::move(mapping), cell.lattice(),
                       tolerance, angle_tolerance};
    }

    auto const prim_lat =
        primitive_lattice(cell, pure, aperiodic_axis, tolerance);
    if (!prim_lat) {
      continue;
    }
    auto trimmed = trim_cell(*prim_lat, cell, aperiodic_axis, tolerance);
    if (!trimmed) {
      continue;
    }
    return Primitive{std::move(trimmed->first), std::move(trimmed->second),
                     cell.lattice(), tolerance, angle_tolerance};
  }
  return leaf::new_error(e_cell_standardization_failed{});
}

std::optional<Matrix3d> primitive_lattice_vectors(
    Cell const &cell, std::vector<Vector3d> const &pure_trans, double symprec) {
  return primitive_lattice(cell, pure_trans, cell.aperiodic_axis(), symprec);
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
    Matrix3d const rot_d = t_mat * op.rotation.cast<double>() * (*t_mat_inv);
    op.rotation = math::round_to_int(rot_d);
    op.translation = t_mat * op.translation;
  }
  return std::make_pair(std::move(*prim), t_mat);
}

Result<Primitive>
find_primitive_with_pure_translations(Cell const &cell,
                                      std::vector<Vector3d> const &pure_trans,
                                      double symprec) {
  std::optional<int> const aperiodic_axis = cell.aperiodic_axis();
  if (pure_trans.size() == 1) {
    auto smallest = smallest_lattice_cell(cell, aperiodic_axis, symprec);
    if (!smallest) {
      return leaf::new_error(e_cell_standardization_failed{});
    }
    std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
    std::iota(mapping.begin(), mapping.end(), 0);
    return Primitive{std::move(*smallest), std::move(mapping), cell.lattice(),
                     symprec, std::nullopt};
  }
  auto const prim_lat =
      primitive_lattice(cell, pure_trans, aperiodic_axis, symprec);
  if (!prim_lat) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  auto trimmed = trim_cell(*prim_lat, cell, aperiodic_axis, symprec);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return Primitive{std::move(trimmed->first), std::move(trimmed->second),
                   cell.lattice(), symprec, std::nullopt};
}

} // namespace spglib::symmetry
