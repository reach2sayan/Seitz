#pragma once

// Bridging helpers between CppCrystal's Eigen value types and the reference
// spglib C API (target Spglib::symspg). Compiled only into the oracle test
// targets, gated by SPGLIB_BUILD_ORACLE_TESTS.
//
// FOOTGUN: spglib's C `double lattice[3][3]` is row-major with the basis
// vectors stored as COLUMNS (a = {L[0][0], L[1][0], L[2][0]}). Our Matrix3d
// stores the basis vectors as columns too and indexes [row][col], so the bridge
// is a plain element-wise copy `L[i][j] = M(i, j)`. Never
// `Eigen::Map<Matrix3d>(&L[0][0])` (column-major map of row-major data silently
// transposes).

#include <spglib/core/cell.hpp>
#include <spglib/core/symmetry_operation.hpp>

extern "C" {
#include <spglib.h>
}

// The standalone spg_get_* functions are deprecated in favour of SpglibDataset,
// but remain the most direct oracle for individual quantities. Silence the
// noise.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif

#include <array>
#include <string>
#include <vector>

namespace spglib::oracle {

inline void to_c_lattice(double out[3][3], Matrix3d const &m) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      out[i][j] = m(i, j);
}

// Owns the C-array storage for a cell so the pointers handed to spglib stay
// valid.
struct CCell {
  double lattice[3][3];
  std::vector<std::array<double, 3>> position;
  std::vector<int> types;

  explicit CCell(Cell const &cell)
      : position(static_cast<std::size_t>(cell.size())),
        types(static_cast<std::size_t>(cell.size())) {
    to_c_lattice(lattice, cell.lattice());
    for (Index i = 0; i < cell.size(); ++i) {
      auto const u = static_cast<std::size_t>(i);
      position[u] = {cell.positions()(i, 0), cell.positions()(i, 1),
                     cell.positions()(i, 2)};
      types[u] = cell.types()[u];
    }
  }

  [[nodiscard]] int num_atom() const { return static_cast<int>(types.size()); }
  [[nodiscard]] double const (*pos() const)[3] {
    return reinterpret_cast<double const(*)[3]>(position.data());
  }
  [[nodiscard]] double (*pos_mut())[3] {
    return reinterpret_cast<double(*)[3]>(position.data());
  }
};

// Symmetry operations as found by the reference spg_get_symmetry, in the input
// cell's basis (rotations integer, translations fractional).
inline SymmetryOperations reference_symmetry(Cell const &cell, double symprec) {
  CCell c(cell);
  int const n = c.num_atom();
  int const mult =
      spg_get_multiplicity(c.lattice, c.pos(), c.types.data(), n, symprec);
  std::vector<int> rot(static_cast<std::size_t>(9 * mult));
  std::vector<double> trans(static_cast<std::size_t>(3 * mult));
  int const found =
      spg_get_symmetry(reinterpret_cast<int(*)[3][3]>(rot.data()),
                       reinterpret_cast<double(*)[3]>(trans.data()), mult,
                       c.lattice, c.pos(), c.types.data(), n, symprec);
  SymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(found));
  for (int s = 0; s < found; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = rot[static_cast<std::size_t>(9 * s + 3 * i + j)];
    Vector3d t(trans[static_cast<std::size_t>(3 * s + 0)],
               trans[static_cast<std::size_t>(3 * s + 1)],
               trans[static_cast<std::size_t>(3 * s + 2)]);
    ops.push_back({r, t});
  }
  return ops;
}

// Reference international space-group number (0 on failure).
inline int reference_spacegroup_number(Cell const &cell, double symprec) {
  CCell c(cell);
  char symbol[11];
  return spg_get_international(symbol, c.lattice, c.pos(), c.types.data(),
                               c.num_atom(), symprec);
}

inline void from_c_lattice(Matrix3d &m, double const in[3][3]) {
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      m(i, j) = in[i][j];
}

// Reference Niggli-reduced lattice (spg_niggli_reduce overwrites in place).
inline Matrix3d reference_niggli(Matrix3d const &lattice, double symprec) {
  double l[3][3];
  to_c_lattice(l, lattice);
  spg_niggli_reduce(l, symprec);
  Matrix3d out;
  from_c_lattice(out, l);
  return out;
}

// Reference Delaunay-reduced lattice (spg_delaunay_reduce overwrites in place).
inline Matrix3d reference_delaunay(Matrix3d const &lattice, double symprec) {
  double l[3][3];
  to_c_lattice(l, lattice);
  spg_delaunay_reduce(l, symprec);
  Matrix3d out;
  from_c_lattice(out, l);
  return out;
}

// Reference primitive cell (spg_find_primitive overwrites in place and returns
// the new atom count).
inline Cell reference_find_primitive(Cell const &cell, double symprec) {
  CCell c(cell);
  int const new_n = spg_find_primitive(c.lattice, c.pos_mut(), c.types.data(),
                                       c.num_atom(), symprec);
  Matrix3d lattice;
  from_c_lattice(lattice, c.lattice);
  Positions pos(new_n, 3);
  Types types(static_cast<std::size_t>(new_n));
  for (int i = 0; i < new_n; ++i) {
    pos.row(i) << c.position[static_cast<std::size_t>(i)][0],
        c.position[static_cast<std::size_t>(i)][1],
        c.position[static_cast<std::size_t>(i)][2];
    types[static_cast<std::size_t>(i)] = c.types[static_cast<std::size_t>(i)];
  }
  return Cell(lattice, pos, types);
}

// Reference database operations for a Hall number
// (spg_get_symmetry_from_database).
inline SymmetryOperations reference_database_operations(int hall_number) {
  int rot[192][3][3];
  double trans[192][3];
  int const n = spg_get_symmetry_from_database(rot, trans, hall_number);
  SymmetryOperations ops;
  ops.reserve(static_cast<std::size_t>(n));
  for (int s = 0; s < n; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = rot[s][i][j];
    ops.push_back({r, Vector3d(trans[s][0], trans[s][1], trans[s][2])});
  }
  return ops;
}

inline SpglibSpacegroupType reference_spacegroup_type(int hall_number) {
  return spg_get_spacegroup_type(hall_number);
}

// Reference space-group determination (spg_get_dataset). number == 0 signals
// failure. Strings are copied out before the dataset is freed.
struct RefDataset {
  int number = 0;
  int hall_number = 0;
  std::string international; // short Hermann-Mauguin
  std::string hall_symbol;
  std::string choice;
  std::string pointgroup;
  int n_operations = 0;
  SymmetryOperations operations;
  int n_std_atoms = 0;
  Matrix3d std_lattice;
  Matrix3d std_rotation_matrix;
  Matrix3d transformation_matrix;
  Vector3d origin_shift;
  // standardized cell positions/types and per-atom Wyckoff data
  Positions std_positions;
  std::vector<int> std_types;
  std::vector<int> wyckoffs;
  std::vector<int> equivalent_atoms;
  std::vector<std::string> site_symmetry_symbols;
};

inline RefDataset reference_dataset(Cell const &cell, double symprec) {
  CCell c(cell);
  SpglibDataset *ds = spg_get_dataset(c.lattice, c.pos_mut(), c.types.data(),
                                      c.num_atom(), symprec);
  RefDataset out;
  if (ds == nullptr) {
    return out;
  }
  out.number = ds->spacegroup_number;
  out.hall_number = ds->hall_number;
  out.international = ds->international_symbol;
  out.hall_symbol = ds->hall_symbol;
  out.choice = ds->choice;
  out.pointgroup = ds->pointgroup_symbol;
  out.n_operations = ds->n_operations;
  for (int s = 0; s < ds->n_operations; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = ds->rotations[s][i][j];
    out.operations.push_back({r, Vector3d(ds->translations[s][0],
                                          ds->translations[s][1],
                                          ds->translations[s][2])});
  }
  from_c_lattice(out.std_lattice, ds->std_lattice);
  from_c_lattice(out.std_rotation_matrix, ds->std_rotation_matrix);
  from_c_lattice(out.transformation_matrix, ds->transformation_matrix);
  out.origin_shift =
      Vector3d(ds->origin_shift[0], ds->origin_shift[1], ds->origin_shift[2]);
  out.n_std_atoms = ds->n_std_atoms;
  out.std_positions = Positions(ds->n_std_atoms, 3);
  out.std_types.resize(static_cast<std::size_t>(ds->n_std_atoms));
  for (int i = 0; i < ds->n_std_atoms; ++i) {
    out.std_positions.row(i) << ds->std_positions[i][0],
        ds->std_positions[i][1], ds->std_positions[i][2];
    out.std_types[static_cast<std::size_t>(i)] = ds->std_types[i];
  }
  out.wyckoffs.assign(ds->wyckoffs, ds->wyckoffs + ds->n_atoms);
  out.equivalent_atoms.assign(ds->equivalent_atoms,
                              ds->equivalent_atoms + ds->n_atoms);
  for (int i = 0; i < ds->n_atoms; ++i) {
    out.site_symmetry_symbols.emplace_back(ds->site_symmetry_symbols[i]);
  }
  spg_free_dataset(ds);
  return out;
}

struct CPointgroup {
  int number;
  std::string symbol;
  Matrix3i transformation;
};

// Reference point group of a set of rotations (spg_get_pointgroup).
inline CPointgroup
reference_pointgroup(std::vector<Matrix3i> const &rotations) {
  std::vector<int> rot(rotations.size() * 9);
  for (std::size_t s = 0; s < rotations.size(); ++s)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        rot[s * 9 + static_cast<std::size_t>(3 * i + j)] = rotations[s](i, j);
  char symbol[6] = {};
  int tmat[3][3];
  int const num = spg_get_pointgroup(symbol, tmat,
                                     reinterpret_cast<int(*)[3][3]>(rot.data()),
                                     static_cast<int>(rotations.size()));
  Matrix3i t;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      t(i, j) = tmat[i][j];
  return {num, std::string(symbol), t};
}

} // namespace spglib::oracle
