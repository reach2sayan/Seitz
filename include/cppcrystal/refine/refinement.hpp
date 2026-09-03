#pragma once

#include <cppcrystal/core/types.hpp>
#include <cppcrystal/spacegroup/spacegroup.hpp>

// Standardization (3D space-group path): idealized conventional lattice, the
// rigid rotation that orients it, and the adjustment of the bravais lattice /
// origin shift to the setting most similar to the idealized one. The
// exact-positions / Wyckoff assembly builds on these.
namespace cppcrystal::refine {

// Idealized conventional lattice built purely from the metric (lengths +
// angles) of the spacegroup's bravais lattice, in the canonical orientation for
// the crystal system.
[[nodiscard]] Matrix3d conventional_lattice(spacegroup::Spacegroup const &sg);

// Rotate the spacegroup's bravais lattice — and correspondingly its origin
// shift — to the proper-rotation setting whose basis vectors are closest
// (Frobenius) to the idealized conventional lattice. Returns the adjusted
// Spacegroup.
[[nodiscard]] spacegroup::Spacegroup
find_similar_bravais_lattice(spacegroup::Spacegroup sg, double symprec);

} // namespace cppcrystal::refine
