#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>

#include <vector>

namespace spglib::symmetry {

struct Primitive {
  Cell cell;                      // primitive cell (Delaunay-reduced lattice)
  std::vector<int> mapping_table; // input atom index -> primitive atom index
};

// Find the primitive cell of `cell`. Port of primitive.c prm_get_primitive /
// the public spg_find_primitive. Errors with e_cell_standardization_failed when
// no primitive cell can be determined.
[[nodiscard]] Result<Primitive>
find_primitive(Cell const &cell, double symprec,
               AngleTolerance angle_tolerance = std::nullopt);

} // namespace spglib::symmetry
