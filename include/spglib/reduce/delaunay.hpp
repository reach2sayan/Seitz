#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/types.hpp>

namespace spglib::reduce {

// Delaunay reduction of a 3D lattice (columns = basis vectors). Returns the
// reduced lattice (right-handed, columns = reduced basis vectors). Port of
// spglib's del_delaunay_reduce (delaunay.c), space-group path only.
//
// Errors with e_delaunay_failed when the cell is degenerate or the change of
// basis is not unimodular.
[[nodiscard]] Result<Matrix3d> delaunay_reduce(Matrix3d const &lattice,
                                               double symprec);

} // namespace spglib::reduce
