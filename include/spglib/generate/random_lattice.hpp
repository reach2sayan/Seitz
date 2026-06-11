#pragma once

#include <spglib/core/types.hpp>

#include <cstdint>

namespace spglib::generate {

// The seven crystal systems, which fix the free lattice parameters. Derived from
// the international space-group number, not the centering (the centering selects
// the Bravais lattice within a system but does not change which of a, b, c,
// alpha, beta, gamma are independent).
enum class CrystalSystem {
  triclinic,
  monoclinic,
  orthorhombic,
  tetragonal,
  trigonal,
  hexagonal,
  cubic,
};

// The crystal system of an international space-group number (1..230), by the
// standard contiguous ranges.
[[nodiscard]] CrystalSystem crystal_system(int spacegroup_number) noexcept;

// A random lattice (columns = Cartesian basis vectors) consistent with the
// crystal system's metric constraints, scaled to `target_volume`. Deterministic
// in `seed`. The unconstrained parameters are sampled from sensible ranges
// (lengths around unity before scaling, angles within the usual crystallographic
// bounds); the absolute scale is then fixed by `target_volume`.
[[nodiscard]] Matrix3d random_lattice(CrystalSystem system, double target_volume,
                                      std::uint64_t seed);

} // namespace spglib::generate
