#pragma once

#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>
#include <seitz/generate/assignments.hpp> // Composition

#include <cstdint>
#include <span>

// The random metrics the generators draw their cells from. Private: a caller
// asks a Generator for a structure, not for a lattice.
namespace seitz::generate {

enum class CrystalSystem {
  triclinic,
  monoclinic,
  orthorhombic,
  tetragonal,
  trigonal,
  hexagonal,
  cubic,
};

[[nodiscard]] constexpr CrystalSystem
crystal_system(int spacegroup_number) noexcept {
  int const n = spacegroup_number;
  if (n <= 2) {
    return CrystalSystem::triclinic;
  } else if (n <= 15) {
    return CrystalSystem::monoclinic;
  } else if (n <= 74) {
    return CrystalSystem::orthorhombic;
  } else if (n <= 142) {
    return CrystalSystem::tetragonal;
  } else if (n <= 167) {
    return CrystalSystem::trigonal;
  } else if (n <= 194) {
    return CrystalSystem::hexagonal;
  }
  return CrystalSystem::cubic;
}

// A random lattice (columns = Cartesian basis vectors) consistent with the
// crystal system's metric constraints, scaled to `target_volume`. Deterministic
// in `seed`. The unconstrained parameters are sampled from sensible ranges
// (lengths around unity before scaling, angles within the usual
// crystallographic bounds); the absolute scale is then fixed by
// `target_volume`.
[[nodiscard]] Matrix3d random_lattice(CrystalSystem system,
                                      double target_volume, std::uint64_t seed);

// Element-aware estimate of the conventional-cell volume (A^3) of a
// composition: each atom contributes (4/3) pi r_cov^3 (data::atomic_volume),
// the sum divided by a representative packing fraction so it lands near real
// cell volumes rather than the cramped sphere-sum. Untabulated types use
// `fallback_volume`.
[[nodiscard]] double estimated_cell_volume(Composition const &composition,
                                           double fallback_volume = 20.0);

// A random layer lattice (columns = Cartesian basis vectors): a, b span the
// periodic plane, c is perpendicular with length `c_length` (aperiodic axis 2).
// The in-plane metric is symmetrized over the 2x2 blocks of `operations`, so it
// is exactly invariant under the in-plane point group whatever the crystal
// system -- a number-range table would mishandle mixed cases (p112 wants an
// oblique gamma where pm11 forces gamma = 90 in the same range). Scaled to
// `target_area`; deterministic in `seed`.
[[nodiscard]] Matrix3d
random_layer_lattice(std::span<SymmetryOperation const> operations,
                     double target_area, double c_length, std::uint64_t seed);

} // namespace seitz::generate
