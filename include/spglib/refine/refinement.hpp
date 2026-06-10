#pragma once

#include <spglib/core/types.hpp>
#include <spglib/spacegroup/spacegroup.hpp>

// Standardization (refinement.c, 3D space-group path): idealized conventional
// lattice, the rigid rotation that orients it, and the adjustment of the
// bravais lattice / origin shift to the setting most similar to the idealized
// one. The exact-positions / Wyckoff assembly builds on these.
namespace spglib::refine {

// Idealized conventional lattice built purely from the metric (lengths +
// angles) of the spacegroup's bravais lattice, in the canonical orientation for
// the crystal system. Port of refinement.c ref_get_conventional_lattice +
// set_{tricli,monocli,ortho,tetra,trigo,rhomb,cubic}.
[[nodiscard]] Matrix3d conventional_lattice(spacegroup::Spacegroup const &sg);

// Rotate the spacegroup's bravais lattice — and correspondingly its origin
// shift — to the proper-rotation setting whose basis vectors are closest
// (Frobenius) to the idealized conventional lattice. Returns the adjusted
// Spacegroup. Port of refinement.c ref_find_similar_bravais_lattice.
[[nodiscard]] spacegroup::Spacegroup
find_similar_bravais_lattice(spacegroup::Spacegroup sg, double symprec);

// Rigid rotation R such that std_lattice = R . bravais_lattice, built from the
// orthonormal frames of the two lattices' first two basis vectors. Port of
// refinement.c measure_rigid_rotation / get_orthonormal_basis.
[[nodiscard]] Matrix3d measure_rigid_rotation(Matrix3d const &bravais_lattice,
                                              Matrix3d const &std_lattice);

} // namespace spglib::refine
