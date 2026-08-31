#include <cppcrystal/generate/random_lattice.hpp>

#include <cppcrystal/data/element_data.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>

namespace cppcrystal::generate {

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

double estimated_cell_volume(std::map<int, int> const &composition,
                             double fallback_volume) {
  // Representative atomic packing fraction: covalent-radius spheres fill roughly
  // this fraction of a real cell, so dividing the sphere-sum by it recovers a
  // realistic cell volume (e.g. rock-salt NaCl lands near its true ~179 A^3).
  constexpr double kPackingFraction = 0.55;

  double const sphere_sum = std::ranges::fold_left(
      composition, 0.0, [&](double sum, auto const &entry) {
        auto const &[type, count] = entry;
        double const v = data::atomic_volume(type).value_or(fallback_volume);
        return sum + static_cast<double>(count) * v;
      });
  return sphere_sum / kPackingFraction;
}

Matrix3d random_layer_lattice(std::span<SymmetryOperation const> operations,
                              double target_area, double c_length,
                              std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> len(1.0, 2.0);
  std::uniform_real_distribution<double> ang(70.0, 110.0);

  // A generic positive-definite in-plane metric, then averaged over the in-plane
  // 2x2 blocks of the operations so it is exactly invariant under the in-plane
  // point group (whatever crystal system that implies).
  double const a0 = len(rng);
  double const b0 = len(rng);
  double const g0 = ang(rng) * kDeg;
  Eigen::Matrix2d g0m;
  g0m << a0 * a0, a0 * b0 * std::cos(g0), a0 * b0 * std::cos(g0), b0 * b0;

  Eigen::Matrix2d metric = Eigen::Matrix2d::Zero();
  for (auto const &op : operations) {
    Eigen::Matrix2d const r = op.rotation.topLeftCorner<2, 2>().cast<double>();
    metric += r.transpose() * g0m * r;
  }
  if (operations.empty()) {
    metric = g0m;
  } else {
    metric /= static_cast<double>(operations.size());
  }

  // Recover a, b, gamma from the symmetrized metric.
  double a = std::sqrt(metric(0, 0));
  double b = std::sqrt(metric(1, 1));
  double const cos_g = metric(0, 1) / (a * b);
  double const sin_g = std::sqrt(std::max(0.0, 1.0 - cos_g * cos_g));

  // Scale the in-plane cell to the requested 2D area (= a b sin gamma).
  double const scale = std::sqrt(target_area / (a * b * sin_g));
  a *= scale;
  b *= scale;

  Matrix3d lattice;
  lattice.col(0) = Vector3d{a, 0.0, 0.0};
  lattice.col(1) = Vector3d{b * cos_g, b * sin_g, 0.0};
  lattice.col(2) = Vector3d{0.0, 0.0, c_length}; // perpendicular, non-periodic
  return lattice;
}

} // namespace cppcrystal::generate
