#include <spglib/generate/distance_check.hpp>

#include <spglib/data/element_data.hpp>
#include <spglib/math/integer_matrix.hpp>

#include <algorithm>
#include <limits>
#include <vector>

namespace spglib::generate {

namespace {

// Covalent radius of an atom type, falling back to tol.fallback_radius for an
// untabulated element so the check stays well-defined. constexpr — a pure lookup
// over the compile-time covalent-radius table.
[[nodiscard]] constexpr double radius_of(int type,
                                         DistanceTolerance const &tol) noexcept {
  return data::covalent_radius(type).value_or(tol.fallback_radius);
}

} // namespace

double minimum_image_distance(Vector3d const &a, Vector3d const &b,
                              Matrix3d const &lattice, bool include_origin,
                              std::optional<int> aperiodic_axis) noexcept {
  // Fold the fractional difference into [-0.5, 0.5]^3, then search the adjacent
  // cells: for skewed lattices the true nearest image can sit one cell over.
  Vector3d const diff = a - b;
  Vector3d base = math::nearest_offset(diff);
  if (aperiodic_axis) {
    base[*aperiodic_axis] = diff[*aperiodic_axis]; // not periodic: keep raw
  }

  double best = std::numeric_limits<double>::infinity();
  for (int n0 = -1; n0 <= 1; ++n0) {
    for (int n1 = -1; n1 <= 1; ++n1) {
      for (int n2 = -1; n2 <= 1; ++n2) {
        Vector3d neighbour(n0, n1, n2);
        if (aperiodic_axis) {
          neighbour[*aperiodic_axis] = 0.0; // no periodic images along it
        }
        Vector3d const off = base + neighbour;
        if (!include_origin && off.squaredNorm() < 1e-20) {
          continue; // the home image of an atom against itself
        }
        double const d2 = (lattice * off).squaredNorm();
        if (d2 < best) {
          best = d2;
        }
      }
    }
  }
  return std::sqrt(best);
}

bool distances_valid(Cell const &cell, DistanceTolerance tol) noexcept {
  Index const n = cell.size();
  std::vector<double> radius(static_cast<std::size_t>(n));
  std::ranges::transform(cell.types(), radius.begin(),
                         [&](int type) { return radius_of(type, tol); });

  for (Index i = 0; i < n; ++i) {
    Vector3d const pi = cell.position(i);
    double const ri = radius[static_cast<std::size_t>(i)];
    // j == i covers the atom against its own periodic images (cell-too-small);
    // j > i covers distinct pairs and their images.
    for (Index j = i; j < n; ++j) {
      double const min_dist =
          tol.scale * (ri + radius[static_cast<std::size_t>(j)]);
      double const dist = minimum_image_distance(
          pi, cell.position(j), cell.lattice(), /*include_origin=*/i != j,
          cell.aperiodic_axis());
      if (dist < min_dist) {
        return false;
      }
    }
  }
  return true;
}

} // namespace spglib::generate
