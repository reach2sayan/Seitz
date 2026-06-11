#include <spglib/generate/random_lattice.hpp>

#include <cmath>
#include <numbers>
#include <random>

namespace spglib::generate {

namespace {

constexpr double kDeg = std::numbers::pi / 180.0;

// Build a lattice (columns = basis vectors) from cell parameters, in the
// standard crystallographic orientation: a along x, b in the xy-plane.
[[nodiscard]] Matrix3d from_parameters(double a, double b, double c,
                                       double alpha, double beta,
                                       double gamma) {
  double const ca = std::cos(alpha * kDeg);
  double const cb = std::cos(beta * kDeg);
  double const cg = std::cos(gamma * kDeg);
  double const sg = std::sin(gamma * kDeg);

  Vector3d const av{a, 0.0, 0.0};
  Vector3d const bv{b * cg, b * sg, 0.0};
  double const cx = c * cb;
  double const cy = c * (ca - cb * cg) / sg;
  double const cz = c *
                    std::sqrt(std::max(0.0, 1.0 - ca * ca - cb * cb - cg * cg +
                                                2.0 * ca * cb * cg)) /
                    sg;
  Vector3d const cv{cx, cy, cz};

  Matrix3d lattice;
  lattice.col(0) = av;
  lattice.col(1) = bv;
  lattice.col(2) = cv;
  return lattice;
}

} // namespace

CrystalSystem crystal_system(int spacegroup_number) noexcept {
  int const n = spacegroup_number;
  if (n <= 2) {
    return CrystalSystem::triclinic;
  }
  if (n <= 15) {
    return CrystalSystem::monoclinic;
  }
  if (n <= 74) {
    return CrystalSystem::orthorhombic;
  }
  if (n <= 142) {
    return CrystalSystem::tetragonal;
  }
  if (n <= 167) {
    return CrystalSystem::trigonal;
  }
  if (n <= 194) {
    return CrystalSystem::hexagonal;
  }
  return CrystalSystem::cubic;
}

Matrix3d random_lattice(CrystalSystem system, double target_volume,
                        std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> len(1.0, 2.0);    // relative lengths
  std::uniform_real_distribution<double> ang(75.0, 105.0); // generic angles

  double a = len(rng);
  double b = len(rng);
  double c = len(rng);
  double alpha = 90.0;
  double beta = 90.0;
  double gamma = 90.0;

  switch (system) {
  case CrystalSystem::triclinic:
    alpha = ang(rng);
    beta = ang(rng);
    gamma = ang(rng);
    break;
  case CrystalSystem::monoclinic:
    beta = ang(rng); // unique axis b
    break;
  case CrystalSystem::orthorhombic:
    break;
  case CrystalSystem::tetragonal:
    b = a;
    break;
  case CrystalSystem::trigonal:
  case CrystalSystem::hexagonal:
    b = a;
    gamma = 120.0;
    break;
  case CrystalSystem::cubic:
    b = a;
    c = a;
    break;
  }

  Matrix3d const raw = from_parameters(a, b, c, alpha, beta, gamma);
  double const vol = std::abs(raw.determinant());
  double const scale = std::cbrt(target_volume / vol);
  return raw * scale;
}

} // namespace spglib::generate
