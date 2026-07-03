#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/types.hpp>

namespace cppcrystal::reduce {

// Niggli reduction of a 3D lattice (columns = basis vectors), via the
// Krivy-Gruber algorithm. Returns the Niggli-reduced lattice. Space-group path
// only.
//
// `eps` is the tolerance applied to the metric parameters (typically symprec).
// Errors with e_niggli_failed if it does not converge.
[[nodiscard]] Result<Matrix3d> niggli_reduce(Matrix3d const &lattice,
                                             double eps);

} // namespace cppcrystal::reduce
