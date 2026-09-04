#pragma once

#include "core/testable.hpp"
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/core/types.hpp>

#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace cppcrystal::symmetry {

struct Primitive {
  Cell cell;                      // primitive cell (Delaunay-reduced lattice)
  std::vector<int> mapping_table; // input atom index -> primitive atom index
  // The original input cell's lattice (columns = basis vectors). Used by the
  // space-group search to prefer a conventional setting whose basis vectors
  // resemble the input.
  Matrix3d orig_lattice{Matrix3d::Identity()};
  // The (possibly tightened) tolerance at which the primitive cell was found,
  // carried into the space-group search that follows.
  Tolerance tolerance{};
};

// Finds the primitive cell of one cell at one tolerance. The family is a
// compile-time parameter: the layer path fixes the third basis vector to the
// aperiodic lattice vector and reduces only the periodic plane, which is an
// `if constexpr` branch rather than a runtime test.
//
// Non-owning: `cell` must outlive the finder.
template <GroupFamily F> class CPPCRYSTAL_TESTABLE PrimitiveFinder {
public:
  PrimitiveFinder(Cell const &cell, Tolerance const &tol) noexcept
      : cell_(cell), tol_(tol) {}

  // The primitive cell, progressively loosening the tolerance across attempts.
  // Errors with e_cell_standardization_failed when none can be determined.
  [[nodiscard]] Result<Primitive> find() const;

  // The primitive cell given the pure translations explicitly rather than
  // recomputing them from the symmetry. `pure_trans` must include the zero
  // translation. Used by the magnetic standardization, where the translations
  // come from the magnetic symmetry.
  [[nodiscard]] Result<Primitive>
  from_pure_translations(std::span<Vector3d const> pure_trans) const;

  // Just the primitive lattice spanned by a set of pure translations;
  // std::nullopt when no lattice of the implied multiplicity exists. The result
  // is Delaunay-reduced.
  [[nodiscard]] std::optional<Lattice>
  lattice_from_pure_translations(std::span<Vector3d const> pure_trans) const;

  // Fold the atoms into the given (smaller) lattice, de-duplicating
  // translationally-equivalent atoms and averaging their positions. Returns the
  // trimmed cell and the input->trimmed atom mapping, or std::nullopt if the
  // atoms do not divide evenly into the lattice.
  [[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
  trim_to(Lattice const &trimmed_lattice) const;

private:
  Cell const &cell_;
  Tolerance tol_;
};

extern template class PrimitiveFinder<GroupFamily::space>;
extern template class PrimitiveFinder<GroupFamily::layer>;

} // namespace cppcrystal::symmetry
