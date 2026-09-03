#include <cppcrystal/refine/standardize.hpp>

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/core/position_index.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/refine/refinement.hpp>
#include <cppcrystal/refine/site_symmetry.hpp>

#include <algorithm>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

// Wyckoff-position assembly (3D path):
//   conventional primitive cell -> exact Wyckoff positions -> expansion across
//   the centering translations into the bravais cell -> per-input-atom Wyckoff
//   / equivalence data (including the supercell broken-symmetry case).
namespace cppcrystal::refine {

using data::operations_from_database;

namespace {

// In the conventional/standardized setting a layer group always has its
// aperiodic axis as c (axis 2); a 3D group has none.
[[nodiscard]] std::optional<int> conventional_aperiodic_axis(int hall_number) {
  return hall_number < 0 ? std::optional<int>(2) : std::nullopt;
}

// Number of pure (identity-rotation) translations among the conventional ops.
[[nodiscard]] int num_pure_translations(SymmetryOperations const &conv_sym) {
  return static_cast<int>(std::ranges::count_if(
      conv_sym, &SymmetryOperation::is_identity_rotation));
}

// Primitive atoms expressed wrt the (idealized) conventional lattice, shifted
// by the origin shift and folded into the cell.
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

// Replicate the exact conventional-primitive atoms across the pure translations
// to build the full bravais cell.
struct Bravais {
  Cell cell;
  std::vector<int> std_mapping_to_primitive; // per bravais atom
};

[[nodiscard]] Bravais expand_in_bravais(Cell const &conv_prim,
                                        Matrix3d const &std_lattice,
                                        SymmetryOperations const &conv_sym,
                                        ExactPositions const &exact,
                                        std::optional<int> aperiodic_axis) {
  auto const total =
      exact.size() * static_cast<std::size_t>(num_pure_translations(conv_sym));

  std::vector<Vector3d> pos;
  Types types;
  std::vector<int> mapping;
  pos.reserve(total);
  types.reserve(total);
  mapping.reserve(total);
  for (auto const &op : conv_sym | std::views::filter(
                            &SymmetryOperation::is_identity_rotation)) {
    for (auto const &[j, atom] : exact | std::views::enumerate) {
      pos.emplace_back(
          math::wrap_to_unit_cell(Vector3d(atom.position + op.translation)));
      types.push_back(conv_prim.type(j));
      mapping.push_back(static_cast<int>(j));
    }
  }
  return {
      Cell(std_lattice, to_positions(pos), std::move(types), aperiodic_axis),
      std::move(mapping)};
}

// first_index_of(values, size)[v] = the first position holding value v (for
// v in [0, size)); nullopt for values that never occur.
[[nodiscard]] std::vector<std::optional<int>>
first_index_of(std::vector<int> const &values, std::size_t size) {
  std::vector<std::optional<int>> first(size);
  for (auto const [i, v] : values | std::views::enumerate) {
    auto &slot = first[static_cast<std::size_t>(v)];
    if (!slot) {
      slot = static_cast<int>(i);
    }
  }
  return first;
}

// First atom (per operation, lowest index first) that an operation maps `i`
// onto among the atoms before it; `i` itself when none is found.
[[nodiscard]] int search_equivalent_atom(int i, Cell const &cell,
                                         PositionIndex const &index,
                                         SymmetryOperations const &operations) {
  for (auto const &op : operations) {
    if (auto const j = index.first_match(op.apply(cell.position(i)),
                                         cell.type(i),
                                         [&](int k) { return k < i; })) {
      return *j;
    }
  }
  return i;
}

// Equivalence by the actual found operations (used when the input cell breaks
// the ideal multiplicity). Each atom links through the first atom sharing its
// primitive atom; the class heads chain through the first symmetry image found
// (a first-match chain, deliberately not a full connected-components closure).
[[nodiscard]] std::vector<int>
equivalent_atoms_broken(Cell const &cell, SymmetryOperations const &operations,
                        std::vector<int> const &mapping_table, double symprec) {
  PositionIndex const index(cell, symprec);
  auto const first_sharing = first_index_of(mapping_table, mapping_table.size());
  std::vector<int> equiv;
  equiv.reserve(mapping_table.size());
  for (auto const [i, prim] : mapping_table | std::views::enumerate) {
    // First atom sharing i's primitive atom; i itself when i is the first.
    int const first = *first_sharing[static_cast<std::size_t>(prim)];
    if (first != i) {
      equiv.push_back(equiv[static_cast<std::size_t>(first)]);
      continue;
    }
    int const found =
        search_equivalent_atom(static_cast<int>(i), cell, index, operations);
    equiv.push_back(found == i ? static_cast<int>(i)
                               : equiv[static_cast<std::size_t>(found)]);
  }
  return equiv;
}

// Representative input-cell atom per atom, derived from the primitive-cell
// equivalence classes.
[[nodiscard]] std::vector<int>
crystallographic_orbits(ExactPositions const &exact,
                        std::vector<int> const &mapping_table) {
  // For each primitive atom, the first input-cell atom in its orbit.
  auto const first_of = first_index_of(mapping_table, exact.size());
  std::vector<int> rep;
  rep.reserve(exact.size());
  std::ranges::transform(exact, std::back_inserter(rep),
                         [&](auto const &atom) {
                           return first_of[static_cast<std::size_t>(
                                               atom.equivalent_atom)]
                               .value_or(0);
                         });

  std::vector<int> orbits;
  orbits.reserve(mapping_table.size());
  for (int const prim : mapping_table) {
    orbits.push_back(rep[static_cast<std::size_t>(prim)]);
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

  Bravais bravais =
      expand_in_bravais(conv_prim, std_lattice, conv_sym, *exact, conv_ap);

  // Per input-cell atom Wyckoff letter + site-symmetry symbol (the first
  // bravais block is in primitive-atom order).
  Standardized out;
  out.bravais = std::move(bravais.cell);
  out.std_mapping_to_primitive = std::move(bravais.std_mapping_to_primitive);
  out.wyckoffs.reserve(mapping_table.size());
  out.site_symmetry_symbols.reserve(mapping_table.size());
  for (int const prim : mapping_table) {
    auto const &atom = (*exact)[static_cast<std::size_t>(prim)];
    out.wyckoffs.push_back(atom.wyckoff);
    out.site_symmetry_symbols.push_back(atom.site_symmetry_symbol);
  }

  out.crystallographic_orbits = crystallographic_orbits(*exact, mapping_table);

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

} // namespace cppcrystal::refine
