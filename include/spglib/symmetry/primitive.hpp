#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>

#include <vector>

namespace spglib::symmetry {

struct Primitive {
  Cell cell;                      // primitive cell (Delaunay-reduced lattice)
  std::vector<int> mapping_table; // input atom index -> primitive atom index
  // The original input cell's lattice (columns = basis vectors). Used by the
  // space-group search to prefer a conventional setting whose basis vectors
  // resemble the input. (primitive.c Primitive::orig_lattice)
  Matrix3d orig_lattice{Matrix3d::Identity()};
  // The (possibly tightened) symprec at which the primitive cell was found, and
  // the angle tolerance carried into the symmetry search.
  double tolerance{0.0};
  AngleTolerance angle_tolerance{std::nullopt};
};

// Find the primitive cell of `cell`. Port of primitive.c prm_get_primitive /
// the public spg_find_primitive. Errors with e_cell_standardization_failed when
// no primitive cell can be determined.
[[nodiscard]] Result<Primitive>
find_primitive(Cell const &cell, double symprec,
               AngleTolerance angle_tolerance = std::nullopt);

} // namespace spglib::symmetry
