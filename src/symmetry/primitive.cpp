#include <spglib/symmetry/primitive.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>
#include <spglib/reduce/delaunay.hpp>
#include <spglib/symmetry/find_symmetry.hpp>

#include <cmath>
#include <optional>
#include <utility>

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

// get_cell_with_smallest_lattice: the multiplicity-one case.
[[nodiscard]] std::optional<Cell> smallest_lattice_cell(Cell const &cell,
                                                        double symprec) {
  auto const min_lat = reduce::delaunay_reduce(cell.lattice(), symprec);
  if (!min_lat)
    return std::nullopt;
  Matrix3d const trans = min_lat->inverse() * cell.lattice();
  Positions pos(cell.size(), 3);
  for (Index i = 0; i < cell.size(); ++i)
    pos.row(i) = math::mod1(Vector3d(trans * cell.position(i))).transpose();
  return Cell(*min_lat, pos, cell.types());
}

// find_primitive_lattice_vectors + cleaning + Delaunay reduce. Picks the first
// triple of {pure translations, unit vectors} spanning the primitive volume.
[[nodiscard]] std::optional<Matrix3d>
primitive_lattice(Cell const &cell, std::vector<Vector3d> const &pure_trans,
                  double symprec) {
  int const multi = static_cast<int>(pure_trans.size());
  std::vector<Vector3d> cand = pure_trans; // includes the zero translation
  cand.push_back(Vector3d::UnitX());
  cand.push_back(Vector3d::UnitY());
  cand.push_back(Vector3d::UnitZ());

  double const init_volume = std::abs(cell.lattice().determinant());
  int const n = static_cast<int>(cand.size());
  auto const at = [&](int i) { return cand[static_cast<std::size_t>(i)]; };

  for (int i = 0; i < n; ++i)
    for (int j = i + 1; j < n; ++j)
      for (int k = j + 1; k < n; ++k) {
        Matrix3d tmp;
        tmp.col(0) = cell.lattice() * at(i);
        tmp.col(1) = cell.lattice() * at(j);
        tmp.col(2) = cell.lattice() * at(k);
        double const volume = std::abs(tmp.determinant());
        if (volume <= symprec || math::nint(init_volume / volume) != multi)
          continue;

        Matrix3d relative;
        relative.col(0) = at(i);
        relative.col(1) = at(j);
        relative.col(2) = at(k);
        // Clean the relative lattice via its exact integer inverse, if
        // possible.
        auto const rel_inv = math::inverse(relative, symprec);
        if (rel_inv) {
          Matrix3i const inv_int = math::round_to_int(*rel_inv);
          if (std::abs(inv_int.determinant()) == multi)
            relative = inv_int.cast<double>().inverse();
        }
        Matrix3d const prim = cell.lattice() * relative;
        auto const reduced = reduce::delaunay_reduce(prim, symprec);
        if (reduced)
          return *reduced;
        return std::nullopt; // matches spglib: first valid triple only
      }
  return std::nullopt;
}

// trim_cell + get_overlap_table: fold atoms into the trimmed lattice, de-dup
// translationally-equivalent ones, average positions.
[[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
trim_cell(Matrix3d const &trimmed_lattice, Cell const &cell, double symprec) {
  int const n = static_cast<int>(cell.size());
  int const ratio = std::abs(
      math::nint(cell.lattice().determinant() / trimmed_lattice.determinant()));
  if (ratio == 0 || n % ratio != 0)
    return std::nullopt;

  Matrix3i const tmat =
      math::round_to_int(Matrix3d(trimmed_lattice.inverse() * cell.lattice()));
  if (std::abs(tmat.determinant()) != ratio)
    return std::nullopt;

  // Atom positions expressed in the trimmed basis, folded into [0, 1).
  Positions pos(n, 3);
  for (int i = 0; i < n; ++i)
    pos.row(i) = math::mod1(Vector3d(tmat.cast<double>() * cell.position(i)))
                     .transpose();

  // Overlap table with tolerance adjustment until each class has `ratio` atoms.
  std::vector<int> overlap(static_cast<std::size_t>(n));
  double tol = symprec;
  bool ok = false;
  for (int attempt = 0; attempt < kTrimNumAttempt && !ok; ++attempt) {
    for (int i = 0; i < n; ++i) {
      overlap[static_cast<std::size_t>(i)] = i;
      for (int j = 0; j < n; ++j)
        if (cell.type(i) == cell.type(j) &&
            overlap[static_cast<std::size_t>(j)] == j &&
            is_overlap(row(pos, i), row(pos, j), trimmed_lattice, tol)) {
          overlap[static_cast<std::size_t>(i)] = j;
          break;
        }
    }
    bool retry = false;
    for (int i = 0; i < n && !retry; ++i) {
      if (overlap[static_cast<std::size_t>(i)] != i)
        continue;
      int count = 0;
      for (int j = 0; j < n; ++j)
        if (overlap[static_cast<std::size_t>(j)] == i)
          ++count;
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
      if (std::abs(pk - pi) > 0.5)
        tpos(j, l) += pi < pk ? pi + 1.0 : pi - 1.0;
      else
        tpos(j, l) += pi;
    }
  }
  int const multi = n / tn;
  for (int i = 0; i < tn; ++i)
    for (int l = 0; l < 3; ++l)
      tpos(i, l) = math::mod1(tpos(i, l) / multi);

  return std::make_pair(Cell(trimmed_lattice, tpos, trimmed_types), mapping);
}

} // namespace

Result<Primitive> find_primitive(Cell const &cell, double symprec,
                                 AngleTolerance /*angle_tolerance*/) {
  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttempt;
       ++attempt, tolerance *= kReduceRate) {
    auto const pure = pure_translations(cell, tolerance);
    if (pure.empty())
      continue;

    if (pure.size() == 1) {
      auto smallest = smallest_lattice_cell(cell, tolerance);
      if (!smallest)
        continue;
      std::vector<int> mapping(static_cast<std::size_t>(cell.size()));
      for (Index i = 0; i < cell.size(); ++i)
        mapping[static_cast<std::size_t>(i)] = static_cast<int>(i);
      return Primitive{std::move(*smallest), std::move(mapping)};
    }

    auto const prim_lat = primitive_lattice(cell, pure, tolerance);
    if (!prim_lat)
      continue;
    auto trimmed = trim_cell(*prim_lat, cell, tolerance);
    if (!trimmed)
      continue;
    return Primitive{std::move(trimmed->first), std::move(trimmed->second)};
  }
  return leaf::new_error(e_cell_standardization_failed{});
}

} // namespace spglib::symmetry
