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

// Reference grid-point index of a grid address (spg_get_dense_grid_point_from_
// address).
inline std::size_t reference_grid_point_from_address(Vector3i const &address,
                                                     Vector3i const &mesh) {
  int const a[3] = {address[0], address[1], address[2]};
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  return spg_get_dense_grid_point_from_address(a, m);
}

struct RefIrMesh {
  std::size_t num_ir = 0;
  std::vector<std::size_t> ir_mapping_table;
  std::vector<Vector3i> grid_address;
};

// Convert a list of integer rotations to the C [][3][3] layout.
inline std::vector<std::array<std::array<int, 3>, 3>>
to_c_rotations(std::vector<Matrix3i> const &rotations) {
  std::vector<std::array<std::array<int, 3>, 3>> out(rotations.size());
  for (std::size_t s = 0; s < rotations.size(); ++s)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        out[s][static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            rotations[s](i, j);
  return out;
}

inline RefIrMesh extract_ir_mesh(std::vector<std::array<int, 3>> const &ga,
                                 std::vector<std::size_t> const &mapping,
                                 std::size_t num_ir) {
  RefIrMesh out;
  out.num_ir = num_ir;
  out.ir_mapping_table = mapping;
  out.grid_address.reserve(ga.size());
  for (auto const &a : ga)
    out.grid_address.emplace_back(a[0], a[1], a[2]);
  return out;
}

// Reference irreducible reciprocal mesh of a structure
// (spg_get_dense_ir_reciprocal_mesh).
inline RefIrMesh reference_ir_reciprocal_mesh(Cell const &cell,
                                              Vector3i const &mesh,
                                              Vector3i const &is_shift,
                                              bool time_reversal,
                                              double symprec) {
  CCell c(cell);
  auto const prod = static_cast<std::size_t>(mesh[0]) *
                    static_cast<std::size_t>(mesh[1]) *
                    static_cast<std::size_t>(mesh[2]);
  std::vector<std::array<int, 3>> ga(prod);
  std::vector<std::size_t> mapping(prod);
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  int const s[3] = {is_shift[0], is_shift[1], is_shift[2]};
  std::size_t const num_ir = spg_get_dense_ir_reciprocal_mesh(
      reinterpret_cast<int(*)[3]>(ga.data()), mapping.data(), m, s,
      time_reversal ? 1 : 0, c.lattice, c.pos(), c.types.data(), c.num_atom(),
      symprec);
  return extract_ir_mesh(ga, mapping, num_ir);
}

// Reference q-stabilized mesh (spg_get_dense_stabilized_reciprocal_mesh).
inline RefIrMesh reference_stabilized_reciprocal_mesh(
    Vector3i const &mesh, Vector3i const &is_shift, bool time_reversal,
    std::vector<Matrix3i> const &rotations,
    std::vector<Vector3d> const &qpoints) {
  auto const prod = static_cast<std::size_t>(mesh[0]) *
                    static_cast<std::size_t>(mesh[1]) *
                    static_cast<std::size_t>(mesh[2]);
  std::vector<std::array<int, 3>> ga(prod);
  std::vector<std::size_t> mapping(prod);
  auto rot = to_c_rotations(rotations);
  std::vector<std::array<double, 3>> q(qpoints.size());
  for (std::size_t k = 0; k < qpoints.size(); ++k)
    q[k] = {qpoints[k][0], qpoints[k][1], qpoints[k][2]};
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  int const s[3] = {is_shift[0], is_shift[1], is_shift[2]};
  std::size_t const num_ir = spg_get_dense_stabilized_reciprocal_mesh(
      reinterpret_cast<int(*)[3]>(ga.data()), mapping.data(), m, s,
      time_reversal ? 1 : 0, static_cast<int>(rotations.size()),
      reinterpret_cast<int(*)[3][3]>(rot.data()),
      static_cast<int>(qpoints.size()),
      reinterpret_cast<double(*)[3]>(q.data()));
  return extract_ir_mesh(ga, mapping, num_ir);
}

// Reference grid points reached from an address by each reciprocal rotation
// (spg_get_dense_grid_points_by_rotations).
inline std::vector<std::size_t> reference_grid_points_by_rotations(
    Vector3i const &address_orig, std::vector<Matrix3i> const &rot_reciprocal,
    Vector3i const &mesh, Vector3i const &is_shift) {
  std::vector<std::size_t> out(rot_reciprocal.size());
  auto rot = to_c_rotations(rot_reciprocal);
  int const a[3] = {address_orig[0], address_orig[1], address_orig[2]};
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  int const s[3] = {is_shift[0], is_shift[1], is_shift[2]};
  spg_get_dense_grid_points_by_rotations(
      out.data(), a, static_cast<int>(rot_reciprocal.size()),
      reinterpret_cast<int(*)[3][3]>(rot.data()), m, s);
  return out;
}

struct RefBzGrid {
  std::size_t num_bzgp = 0;
  std::vector<Vector3i> bz_grid_address;     // size num_bzgp
  std::vector<std::size_t> bz_map;           // size 8*prod(mesh), num_bzmesh = unmapped
  std::size_t num_bzmesh = 0;                // the "unmapped" sentinel value
};

// Reference BZ relocation (spg_relocate_dense_BZ_grid_address). bz_grid_address
// is allocated to the worst case (8*prod(mesh)) and trimmed to num_bzgp.
inline RefBzGrid reference_relocate_BZ(std::vector<Vector3i> const &grid_address,
                                       Vector3i const &mesh,
                                       Matrix3d const &rec_lattice,
                                       Vector3i const &is_shift) {
  auto const num_bzmesh = static_cast<std::size_t>(2 * mesh[0]) *
                          static_cast<std::size_t>(2 * mesh[1]) *
                          static_cast<std::size_t>(2 * mesh[2]);
  std::vector<std::array<int, 3>> ga(grid_address.size());
  for (std::size_t i = 0; i < grid_address.size(); ++i)
    ga[i] = {grid_address[i][0], grid_address[i][1], grid_address[i][2]};
  std::vector<std::array<int, 3>> bz_ga(num_bzmesh);
  std::vector<std::size_t> bz_map(num_bzmesh);
  double rl[3][3];
  to_c_lattice(rl, rec_lattice);
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  int const s[3] = {is_shift[0], is_shift[1], is_shift[2]};
  std::size_t const num_bzgp = spg_relocate_dense_BZ_grid_address(
      reinterpret_cast<int(*)[3]>(bz_ga.data()), bz_map.data(),
      reinterpret_cast<int(*)[3]>(ga.data()), m, rl, s);

  RefBzGrid out;
  out.num_bzgp = num_bzgp;
  out.num_bzmesh = num_bzmesh;
  out.bz_map = bz_map;
  out.bz_grid_address.reserve(num_bzgp);
  for (std::size_t i = 0; i < num_bzgp; ++i)
    out.bz_grid_address.emplace_back(bz_ga[i][0], bz_ga[i][1], bz_ga[i][2]);
  return out;
}

// Reference BZ grid points reached by each reciprocal rotation
// (spg_get_dense_BZ_grid_points_by_rotations).
inline std::vector<std::size_t> reference_BZ_grid_points_by_rotations(
    Vector3i const &address_orig, std::vector<Matrix3i> const &rot_reciprocal,
    Vector3i const &mesh, Vector3i const &is_shift,
    std::vector<std::size_t> const &bz_map) {
  std::vector<std::size_t> out(rot_reciprocal.size());
  auto rot = to_c_rotations(rot_reciprocal);
  int const a[3] = {address_orig[0], address_orig[1], address_orig[2]};
  int const m[3] = {mesh[0], mesh[1], mesh[2]};
  int const s[3] = {is_shift[0], is_shift[1], is_shift[2]};
  spg_get_dense_BZ_grid_points_by_rotations(
      out.data(), a, static_cast<int>(rot_reciprocal.size()),
      reinterpret_cast<int(*)[3][3]>(rot.data()), m, s, bz_map.data());
  return out;
}

} // namespace spglib::oracle
