#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/symmetry_operation.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>

#include <optional>
#include <utility>
#include <vector>

namespace spglib::symmetry {

struct Primitive {
  Cell cell;                      // primitive cell (Delaunay-reduced lattice)
  std::vector<int> mapping_table; // input atom index -> primitive atom index
  // The original input cell's lattice (columns = basis vectors). Used by the
  // space-group search to prefer a conventional setting whose basis vectors
  // resemble the input.
  Matrix3d orig_lattice{Matrix3d::Identity()};
  // The (possibly tightened) symprec at which the primitive cell was found, and
  // the angle tolerance carried into the symmetry search.
  double tolerance{0.0};
  AngleTolerance angle_tolerance{std::nullopt};
};

// Find the primitive cell of `cell`. Errors with e_cell_standardization_failed
// when no primitive cell can be determined.
[[nodiscard]] Result<Primitive>
find_primitive(Cell const &cell, double symprec,
               AngleTolerance angle_tolerance = std::nullopt);

// The primitive lattice (columns = basis vectors) spanned by a set of pure
// translations of `cell` (which must include the zero translation). std::nullopt
// when no primitive lattice of the implied multiplicity (= pure_trans.size())
// exists. The result is Delaunay-reduced. Used by the magnetic-symmetry search
// to rebuild the primitive cell from the magnetic pure translations.
[[nodiscard]] std::optional<Matrix3d>
primitive_lattice_vectors(Cell const &cell,
                          std::vector<Vector3d> const &pure_trans,
                          double symprec);

// The symmetry operations of the primitive cell implied by a set of (typically
// conventional) symmetry operations, together with the transformation `t_mat`
// from the primitive to the conventional setting:
//   (a_p, b_p, c_p) . t_mat = (a_c, b_c, c_c).
// The pure translations of `operations` (rotation == identity) define the
// primitive lattice in "translation space"; the distinct rotations, transformed
// to that setting (R' = T R T^-1, t' = T t), are the primitive operations.
// std::nullopt if the operation set is inconsistent (e.g. its size is not a
// multiple of the number of pure translations, or the translations do not span
// a single primitive point).
[[nodiscard]] std::optional<std::pair<SymmetryOperations, Matrix3d>>
primitive_symmetry(SymmetryOperations const &operations, double symprec);

// The primitive cell of `cell` given its pure translations explicitly (rather
// than recomputing them from the symmetry). `pure_trans` must include the zero
// translation. Used by the magnetic standardization, where the pure
// translations come from the magnetic symmetry. Errors with
// e_cell_standardization_failed.
[[nodiscard]] Result<Primitive>
find_primitive_with_pure_translations(Cell const &cell,
                                      std::vector<Vector3d> const &pure_trans,
                                      double symprec);

// Fold the atoms of `cell` into the given (smaller) `trimmed_lattice`,
// de-duplicating translationally-equivalent atoms and averaging their
// positions. Returns the trimmed cell together with the input->trimmed atom
// mapping, or std::nullopt if the atoms do not divide evenly into the lattice.
// The aperiodic axis (for layer cells) is read from `cell`. Exposed for cell
// standardization (transform-to-primitive).
[[nodiscard]] std::optional<std::pair<Cell, std::vector<int>>>
trim_to_lattice(Matrix3d const &trimmed_lattice, Cell const &cell,
                double symprec);

} // namespace spglib::symmetry
