#pragma once

#include "core/testable.hpp"
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include "refine/standardize.hpp"
#include "spacegroup/spacegroup.hpp"
#include "symmetry/primitive.hpp"

#include <optional>
#include <utility>

// Standardization: the idealized conventional lattice, the rigid rotation that
// orients it, the exact operations recovered in the input cell, and the
// standardized cell with its per-atom Wyckoff data.
namespace cppcrystal::refine {

// Idealized conventional lattice built purely from the metric (lengths +
// angles) of the matched bravais lattice, in the canonical orientation for the
// crystal system. A pure function of the matched group: the magnetic path needs
// it without a cell to refine.
[[nodiscard]] CPPCRYSTAL_TESTABLE Lattice
conventional_lattice(SpacegroupMatch const &sg);

// Rotate the bravais lattice — and correspondingly the origin shift — to the
// proper-rotation setting whose basis vectors are closest (Frobenius) to the
// idealized conventional lattice.
[[nodiscard]] CPPCRYSTAL_TESTABLE SpacegroupMatch
find_similar_bravais_lattice(SpacegroupMatch sg,
                                                           double symprec);

// Turns a matched space group plus the cell it was matched from into the
// standardized result. The family is a compile-time parameter: only the layer
// path leaves the conventional c axis aperiodic.
//
// Non-owning: `primitive` and `cell` must outlive the refinement.
template <GroupFamily F> class CPPCRYSTAL_TESTABLE Refinement {
public:
  Refinement(SpacegroupMatch matched, Cell const &primitive, Cell const &cell,
             Tolerance const &tol)
      : matched_(std::move(matched)), primitive_(primitive), cell_(cell),
        tol_(tol) {}

  // The idealized conventional lattice of the matched group.
  [[nodiscard]] Lattice conventional_lattice() const {
    return refine::conventional_lattice(matched_);
  }

  // The same refinement with the bravais lattice and origin shift rotated to
  // the setting closest to the idealized conventional one.
  [[nodiscard]] Refinement similar_bravais() && {
    return Refinement{
        find_similar_bravais_lattice(std::move(matched_), tol_.symprec),
        primitive_, cell_, tol_};
  }

  [[nodiscard]] SpacegroupMatch const &matched() const noexcept {
    return matched_;
  }

  // Exact, database-derived space-group operations of the input cell: the
  // conventional database operations of the Hall setting, shifted by the origin
  // shift, transformed to the primitive setting, then recovered in the input
  // cell. std::nullopt when the recovered pure-translation count is
  // inconsistent with the input cell's multiplicity (the caller retries at a
  // tighter tolerance). Requires similar_bravais() to have been applied.
  [[nodiscard]] std::optional<Operations> operations() const;

  // Bring `cell` into the primitive setting implied by `transformation` (cell
  // -> conventional) and the matched centering (conventional -> primitive),
  // folding the atoms into the primitive lattice.
  [[nodiscard]] Result<Cell> to_primitive(Cell const &cell,
                                          Matrix3d const &transformation) const;

  // Expand a primitive cell into the centered conventional cell: one copy of
  // each atom per centering translation, folded back into the lattice.
  [[nodiscard]] Result<Cell> from_primitive(Cell const &primitive) const;

  // The standardized cell and the per-atom Wyckoff data. `cell_operations` are
  // the input cell's own operations (used only for the supercell
  // broken-symmetry check); `mapping_table[i]` maps input atom i to its
  // primitive atom. Requires similar_bravais() to have been applied.
  [[nodiscard]] std::optional<Standardized>
  standardize(Operations const &cell_operations,
              std::vector<int> const &mapping_table) const;

private:
  SpacegroupMatch matched_;
  Cell const &primitive_;
  Cell const &cell_;
  Tolerance tol_;
};

// The two members are defined in different translation units, so they are
// instantiated one at a time rather than by whole-class instantiation.
extern template Result<Cell>
Refinement<GroupFamily::space>::to_primitive(Cell const &,
                                             Matrix3d const &) const;
extern template Result<Cell>
Refinement<GroupFamily::layer>::to_primitive(Cell const &,
                                             Matrix3d const &) const;
extern template Result<Cell>
Refinement<GroupFamily::space>::from_primitive(Cell const &) const;
extern template Result<Cell>
Refinement<GroupFamily::layer>::from_primitive(Cell const &) const;
extern template std::optional<Operations>
Refinement<GroupFamily::space>::operations() const;
extern template std::optional<Operations>
Refinement<GroupFamily::layer>::operations() const;
extern template std::optional<Standardized>
Refinement<GroupFamily::space>::standardize(Operations const &,
                                            std::vector<int> const &) const;
extern template std::optional<Standardized>
Refinement<GroupFamily::layer>::standardize(Operations const &,
                                            std::vector<int> const &) const;

} // namespace cppcrystal::refine
