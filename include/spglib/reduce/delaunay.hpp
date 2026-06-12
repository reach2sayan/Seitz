#pragma once

#include <spglib/core/error.hpp>
#include <spglib/core/types.hpp>

namespace spglib::reduce {

// Delaunay reduction of a 3D lattice (columns = basis vectors). Returns the
// reduced lattice (right-handed, columns = reduced basis vectors). Space-group
// path only.
//
// Errors with e_delaunay_failed when the cell is degenerate or the change of
// basis is not unimodular.
[[nodiscard]] Result<Matrix3d> delaunay_reduce(Matrix3d const &lattice,
                                               double symprec);

// 2D Delaunay reduction within the plane spanned by the two axes other than
// `unique_axis`; the unique axis column is left untouched. Space-group path
// only (lattice rank 2). Returns the reduced lattice (columns = basis vectors),
// made right-handed by flipping the unique axis if needed. The extra
// `unique_axis` argument distinguishes this from the 3D overload above.
//
// Errors with e_delaunay_failed when the reduced cell is degenerate.
[[nodiscard]] Result<Matrix3d>
delaunay_reduce(Matrix3d const &lattice, int unique_axis, double symprec);

} // namespace spglib::reduce
