#include <spglib/refine/standardize.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/refine/refinement.hpp>
#include <spglib/refine/site_symmetry.hpp>

#include <algorithm>
#include <optional>

// Port of the get_Wyckoff_positions assembly (refinement.c, 3D path):
//   conventional primitive cell -> exact Wyckoff positions -> expansion across
//   the centering translations into the bravais cell -> per-input-atom Wyckoff
//   / equivalence data (including the supercell broken-symmetry case).
namespace spglib::refine {

using data::operations_from_database;

namespace {

// In the conventional/standardized setting a layer group always has its
// aperiodic axis as c (axis 2); a 3D group has none. Port of refinement.c's
// `conv_prim->aperiodic_axis = hall_number > 0 ? -1 : 2`.
[[nodiscard]] std::optional<int> conventional_aperiodic_axis(int hall_number) {
  return hall_number < 0 ? std::optional<int>(2) : std::nullopt;
}

// Number of pure (identity-rotation) translations among the conventional ops.
[[nodiscard]] int num_pure_translations(SymmetryOperations const &conv_sym) {
  return static_cast<int>(std::ranges::count_if(conv_sym, [](auto const &op) {
    return op.rotation == Matrix3i::Identity();
  }));
}

// get_conventional_primitive: primitive atoms expressed wrt the (idealized)
// conventional lattice, shifted by the origin shift and folded into the cell.
[[nodiscard]] Cell conventional_primitive(spacegroup::Spacegroup const &sg,
                                          Cell const &primitive,
                                          Matrix3d const &std_lattice,
                                          std::optional<int> aperiodic_axis) {
  Matrix3d const trans_mat = sg.bravais_lattice.inverse() * primitive.lattice();
  Positions pos(primitive.size(), 3);
  for (Index i = 0; i < primitive.size(); ++i) {
    Vector3d const p = trans_mat * primitive.position(i) + sg.origin_shift;
    // Positions are stored folded into [0, 1) on every axis — including the
    // aperiodic one: the origin shift has aligned the layer onto the database
    // convention (symmetry plane at c = 0), so folding c resolves the sign
    // ambiguity of the shift. (The aperiodic axis is only left un-folded in the
    // overlap *distance*, via is_overlap, not in stored coordinates.)
    pos.row(i) = math::wrap_to_unit_cell(p).transpose();
  }
  return Cell(std_lattice, std::move(pos), primitive.types(), aperiodic_axis);
}

// expand_positions_in_bravais: replicate the exact conventional-primitive atoms
// across the pure translations to build the full bravais cell.
struct Bravais {
  Cell cell;
  std::vector<int> wyckoffs;
  std::vector<std::string> site_symmetry_symbols;
  std::vector<int> equivalent_atoms;
  std::vector<int> std_mapping_to_primitive;
};

[[nodiscard]] Bravais expand_in_bravais(Cell const &conv_prim,
                                        Matrix3d const &std_lattice,
                                        SymmetryOperations const &conv_sym,
                                        ExactPositions const &exact,
                                        std::optional<int> aperiodic_axis) {
  auto const n = static_cast<std::size_t>(conv_prim.size());
  int const multi = num_pure_translations(conv_sym);
  auto const total = static_cast<Index>(n * static_cast<std::size_t>(multi));

  Positions pos(total, 3);
  Types types(static_cast<std::size_t>(total));
  Bravais b;
  b.wyckoffs.resize(static_cast<std::size_t>(total));
  b.site_symmetry_symbols.resize(static_cast<std::size_t>(total));
  b.equivalent_atoms.resize(static_cast<std::size_t>(total));
  b.std_mapping_to_primitive.resize(static_cast<std::size_t>(total));

  Index idx = 0;
  for (auto const &op : conv_sym) {
    if (op.rotation != Matrix3i::Identity()) {
      continue;
    }
    for (std::size_t j = 0; j < n; ++j) {
      Vector3d const p =
          math::wrap_to_unit_cell(Vector3d(exact.positions[j] + op.translation));
      auto const u = static_cast<std::size_t>(idx);
      pos.row(idx) = p.transpose();
      types[u] = conv_prim.type(static_cast<Index>(j));
      b.wyckoffs[u] = exact.wyckoffs[j];
      b.site_symmetry_symbols[u] = exact.site_symmetry_symbols[j];
      b.equivalent_atoms[u] = exact.equivalent_atoms[j];
      b.std_mapping_to_primitive[u] = static_cast<int>(j);
      ++idx;
    }
  }
  b.cell = Cell(std_lattice, std::move(pos), std::move(types), aperiodic_axis);
  return b;
}

// search_equivalent_atom: lowest-indexed atom that an operation maps `i` onto.
[[nodiscard]] int search_equivalent_atom(int i, Cell const &cell,
                                         SymmetryOperations const &operations,
                                         double symprec) {
  std::optional<int> const aperiodic_axis = cell.aperiodic_axis();
  for (auto const &op : operations) {
    Vector3d const pos = op.rotation.cast<double>() * cell.position(i) +
                         op.translation;
    for (int j = 0; j < i; ++j) {
      if (is_overlap_same_type(cell.position(j), pos, cell.type(j),
                               cell.type(i), cell.lattice(), symprec,
                               aperiodic_axis)) {
        return j;
      }
    }
  }
  return i;
}

// set_equivalent_atoms_broken_symmetry: equivalence by the actual found
// operations (used when the input cell breaks the ideal multiplicity).
[[nodiscard]] std::vector<int>
equivalent_atoms_broken(Cell const &cell, SymmetryOperations const &operations,
                        std::vector<int> const &mapping_table, double symprec) {
  int const n = static_cast<int>(cell.size());
  std::vector<int> equiv(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto const ui = static_cast<std::size_t>(i);
    equiv[ui] = i;
    for (int j = 0; j < n; ++j) {
      if (mapping_table[ui] == mapping_table[static_cast<std::size_t>(j)]) {
        equiv[ui] =
            (i == j)
                ? equiv[static_cast<std::size_t>(
                      search_equivalent_atom(i, cell, operations, symprec))]
                : equiv[static_cast<std::size_t>(j)];
        break;
      }
    }
  }
  return equiv;
}

// set_crystallographic_orbits: representative input-cell atom per atom, derived
// from the primitive-cell equivalence classes.
[[nodiscard]] std::vector<int>
crystallographic_orbits(Cell const &primitive, Cell const &cell,
                        std::vector<int> const &equiv_prim,
                        std::vector<int> const &mapping_table) {
  int const np = static_cast<int>(primitive.size());
  int const nc = static_cast<int>(cell.size());

  // For each primitive atom, the first input-cell atom in its orbit.
  std::vector<int> rep(static_cast<std::size_t>(np));
  for (int i = 0; i < np; ++i) {
    for (int j = 0; j < nc; ++j) {
      if (mapping_table[static_cast<std::size_t>(j)] ==
          equiv_prim[static_cast<std::size_t>(i)]) {
        rep[static_cast<std::size_t>(i)] = j;
        break;
      }
    }
  }

  std::vector<int> orbits(static_cast<std::size_t>(nc));
  for (int i = 0; i < nc; ++i) {
    orbits[static_cast<std::size_t>(i)] =
        rep[static_cast<std::size_t>(mapping_table[static_cast<std::size_t>(i)])];
  }
  return orbits;
}

} // namespace

std::optional<Standardized>
get_wyckoff_positions(spacegroup::Spacegroup const &sg, Cell const &primitive,
                      Cell const &cell,
                      SymmetryOperations const &cell_operations,
                      std::vector<int> const &mapping_table, double symprec) {
  int const hall = sg.type.hall_number;
  SymmetryOperations const conv_sym = operations_from_database(hall);
  int const multi = num_pure_translations(conv_sym);
  std::optional<int> const conv_ap = conventional_aperiodic_axis(hall);

  Matrix3d const std_lattice = conventional_lattice(sg);
  Cell const conv_prim =
      conventional_primitive(sg, primitive, std_lattice, conv_ap);

  auto const exact =
      exact_positions(conv_prim, conv_sym, multi, hall, symprec);
  if (!exact) {
    return std::nullopt;
  }

  Bravais const bravais =
      expand_in_bravais(conv_prim, std_lattice, conv_sym, *exact, conv_ap);

  // Per input-cell atom Wyckoff letter + site-symmetry symbol (the first
  // bravais block is in primitive-atom order).
  auto const nc = static_cast<std::size_t>(cell.size());
  Standardized out;
  out.bravais = bravais.cell;
  out.std_mapping_to_primitive = bravais.std_mapping_to_primitive;
  out.wyckoffs.resize(nc);
  out.site_symmetry_symbols.resize(nc);
  for (std::size_t i = 0; i < nc; ++i) {
    int const prim = mapping_table[i];
    out.wyckoffs[i] = exact->wyckoffs[static_cast<std::size_t>(prim)];
    out.site_symmetry_symbols[i] =
        exact->site_symmetry_symbols[static_cast<std::size_t>(prim)];
  }

  out.crystallographic_orbits = crystallographic_orbits(
      primitive, cell, exact->equivalent_atoms, mapping_table);

  // Equivalent atoms: the crystallographic orbits unless the input cell breaks
  // the ideal site multiplicity (supercell), in which case use the found ops.
  int const num_prim_sym = static_cast<int>(conv_sym.size()) / multi;
  if (static_cast<int>(cell.size()) * num_prim_sym !=
      static_cast<int>(cell_operations.size()) *
          static_cast<int>(primitive.size())) {
    out.equivalent_atoms =
        equivalent_atoms_broken(cell, cell_operations, mapping_table, symprec);
  } else {
    out.equivalent_atoms = out.crystallographic_orbits;
  }

  return out;
}

} // namespace spglib::refine
