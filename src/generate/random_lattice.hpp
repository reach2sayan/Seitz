#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <cstdint>
#include <map>
#include <span>

// The random metrics the generators draw their cells from. Private: a caller
// asks a Generator for a structure, not for a lattice.
namespace cppcrystal::generate {

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
// (lengths around unity before scaling, angles within the usual crystallographic
// bounds); the absolute scale is then fixed by `target_volume`.
[[nodiscard]] Matrix3d random_lattice(CrystalSystem system, double target_volume,
                                      std::uint64_t seed);

// An element-aware estimate of the conventional-cell volume (cubic angstrom)
// for a composition (atom type -> count in the conventional cell). Each atom
// contributes the volume of a sphere of its covalent radius
// (data::atomic_volume); the sum is divided by a representative packing fraction
// so the estimate lands near real cell volumes rather than the cramped
// sphere-sum. Replaces the old flat 20 A^3/atom constant. Untabulated types use
// `fallback_volume` cubic angstrom.
[[nodiscard]] double estimated_cell_volume(std::map<int, int> const &composition,
                                           double fallback_volume = 20.0);

// A random layer lattice (columns = Cartesian basis vectors): a, b span the
// periodic plane, c is perpendicular with length `c_length` (the non-periodic
// stacking direction, aperiodic axis 2). The in-plane metric is symmetrized over
// the in-plane (2x2) blocks of `operations`, so it is exactly invariant under
// the layer group's in-plane point group whatever its crystal system — no
// number-range table, which would mis-handle mixed cases (e.g. p112 wants an
// oblique gamma while pm11 forces gamma = 90 in the same range). Scaled so the
// 2D cell area equals `target_area`. Deterministic in `seed`.
[[nodiscard]] Matrix3d
random_layer_lattice(std::span<SymmetryOperation const> operations,
                     double target_area, double c_length, std::uint64_t seed);

} // namespace cppcrystal::generate
