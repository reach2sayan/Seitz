#pragma once

#include "core/testable.hpp"
#include <seitz/core/cell.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>

#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace seitz::symmetry {

struct Primitive {
  Cell cell;                      // primitive cell (Delaunay-reduced lattice)
  std::vector<int> mapping_table; // input atom index -> primitive atom index
  // The original input cell's lattice (columns = basis vectors). Used by the
  // space-group search to prefer a conventional setting whose basis vectors
  // resemble the input.
  Matrix3d orig_lattice{Matrix3d::Identity()};
  Tolerance tolerance{};
};

// Finds the primitive cell of one cell at one tolerance.
// the layer path fixes the third basis vector to the
// aperiodic lattice vector and reduces only the periodic plane, which is an
// `if constexpr` branch rather than a runtime test.
template <GroupFamily F> class SEITZ_TESTABLE PrimitiveFinder {
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
  // translationally-equivalent atoms and averaging their positions.
  [[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
  trim_to(Lattice const &trimmed_lattice) const;

private:
  Cell const &cell_;
  Tolerance tol_;
};

extern template class PrimitiveFinder<GroupFamily::space>;
extern template class PrimitiveFinder<GroupFamily::layer>;

} // namespace seitz::symmetry
