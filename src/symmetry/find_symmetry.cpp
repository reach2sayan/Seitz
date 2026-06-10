#include <spglib/symmetry/find_symmetry.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>
#include <spglib/reduce/delaunay.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <ranges>

// Port of the symmetry-operation search in symmetry.c (3D space-group path):
//   1. Delaunay-reduce the lattice and enumerate the lattice point group by
//      trying every unimodular integer matrix built from 26 candidate axes that
//      preserves the metric (get_lattice_symmetry).
//   2. Transform those rotations into the input cell's basis
//      (transform_pointsymmetry).
//   3. For each rotation, find the translations that map the cell onto itself
//      (get_translation / get_space_group_operations) via the OverlapChecker.
namespace spglib::symmetry {

namespace {

constexpr int kNumAttempt = 100; // symmetry.c NUM_ATTEMPT
constexpr double kAngleReduceRate = 0.95;
constexpr double kSinDtheta2Cutoff = 1e-12;

// symmetry.c relative_axes: 26 candidate lattice vectors.
constexpr int kRelativeAxes[26][3] = {
    {1, 0, 0},   {0, 1, 0},   {0, 0, 1},  {-1, 0, 0},  {0, -1, 0},
    {0, 0, -1},  {0, 1, 1},   {1, 0, 1},  {1, 1, 0},   {0, -1, -1},
    {-1, 0, -1}, {-1, -1, 0}, {0, 1, -1}, {-1, 0, 1},  {1, -1, 0},
    {0, -1, 1},  {1, 0, -1},  {-1, 1, 0}, {1, 1, 1},   {-1, -1, -1},
    {-1, 1, 1},  {1, -1, 1},  {1, 1, -1}, {1, -1, -1}, {-1, 1, -1},
    {-1, -1, 1},
};

[[nodiscard]] Vector3i axis(int const (&a)[3]) { return {a[0], a[1], a[2]}; }

// symmetry.c is_identity_metric, angle_tolerance < 0 path uses the sin-based
// criterion; a positive angle_tolerance compares angles in degrees.
[[nodiscard]] bool is_identity_metric(Matrix3d const &rotated,
                                      Matrix3d const &orig, double symprec,
                                      AngleTolerance angle_tolerance) {
  const auto length_orig = orig.diagonal().array().sqrt();
  const auto length_rot = rotated.diagonal().array().sqrt();

  auto ok_to_continue =
      ((length_orig - length_rot).abs().maxCoeff() <= symprec);
  if (!ok_to_continue) {
    return false;
  }

  constexpr int elem_sets[3][2] = {{0, 1}, {0, 2}, {1, 2}};
  for (auto const &set : elem_sets) {
    int const j = set[0];
    int const k = set[1];
    if (angle_tolerance) {
      auto angle = [](Matrix3d const &m, int a, int b) {
        return std::acos(m(a, b) / std::sqrt(m(a, a)) / std::sqrt(m(b, b))) /
               std::numbers::pi * 180.0;
      };
      if (std::fabs(angle(orig, j, k) - angle(rotated, j, k)) >
          *angle_tolerance) {
        return false;
      }
    } else {
      double const cos1 = orig(j, k) / length_orig[j] / length_orig[k];
      double const cos2 = rotated(j, k) / length_rot[j] / length_rot[k];
      double const x =
          cos1 * cos2 + std::sqrt(1 - cos1 * cos1) * std::sqrt(1 - cos2 * cos2);
      double const sin_dtheta2 = 1 - x * x;
      double const length_ave2 = ((length_orig[j] + length_rot[j]) *
                                  (length_orig[k] + length_rot[k])) /
                                 4;
      if (sin_dtheta2 > kSinDtheta2Cutoff &&
          sin_dtheta2 * length_ave2 > symprec * symprec) {
        return false;
      }
    }
  }
  return true;
}

// Collect lattice symmetries in the Delaunay basis. std::nullopt signals an
// overflow (> 48 operations), which asks the caller to tighten the tolerance.
[[nodiscard]] std::optional<PointSymmetry>
collect_metric_symmetries(Matrix3d const &min_lattice,
                          Matrix3d const &metric_orig, double symprec,
                          AngleTolerance angle_tolerance) {
  PointSymmetry found;
  for (auto const &ai : kRelativeAxes) {
    for (auto const &aj : kRelativeAxes) {
      for (auto const &ak : kRelativeAxes) {
        Matrix3i axes;
        axes << axis(ai), axis(aj), axis(ak);
        if (int const det = axes.determinant(); det != 1 && det != -1) {
          continue;
        }
        Matrix3d const lattice = min_lattice * axes.cast<double>();
        Matrix3d const metric = math::metric_tensor(lattice);
        if (is_identity_metric(metric, metric_orig, symprec, angle_tolerance)) {
          if (found.size() >= 48) {
            return std::nullopt;
          }
          found.push_back(axes);
        }
      }
    }
  }
  return found;
}

// transform_pointsymmetry: bring rotations from the Delaunay basis
// (original_lat) into the cell basis (new_lat) by the similarity transform b^-1
// . R . b with b = original_lat^-1 . new_lat.
[[nodiscard]] Result<PointSymmetry>
transform_pointsymmetry(PointSymmetry const &orig, Matrix3d const &new_lat,
                        Matrix3d const &original_lat) {
  Matrix3d const trans_mat = original_lat.inverse() * new_lat;
  Matrix3d const trans_inv = trans_mat.inverse();
  double const tol = std::abs(trans_mat.determinant()) / 10.0;

  PointSymmetry out;
  for (Matrix3i const &rot : orig) {
    Matrix3d const similar = trans_inv * rot.cast<double>() * trans_mat;
    if (!math::is_int_matrix(similar, tol)) {
      continue;
    }
    Matrix3i const r = math::round_to_int(similar);
    if (std::abs(r.determinant()) != 1) {
      return leaf::new_error(e_symmetry_operation_search_failed{});
    }
    out.push_back(std::move(r));
  }
  return out;
}

// First-occurrence index of the least-frequent atom type, or nullopt for an
// empty cell (symmetry.c get_index_with_least_atoms).
[[nodiscard]] std::optional<int> index_with_least_atoms(Cell const &cell) {
  int const n = static_cast<int>(cell.size());
  if (n == 0)
    return std::nullopt;

  // How many atoms share each type.
  std::map<int, int> frequency;
  for (int i = 0; i < n; ++i)
    ++frequency[cell.type(i)];

  // Return the first atom whose type is the rarest. Scanning in index order
  // reproduces spglib's tie-break: among equally-rare types, smallest index.
  int const min_frequency = std::ranges::min(frequency | std::views::values);
  for (int i = 0; i < n; ++i)
    if (frequency.at(cell.type(i)) == min_frequency)
      return i;

  return std::nullopt; // unreachable: min_frequency came from the cell itself
}

// Translations t such that x -> rot . x + t maps the cell onto itself
// (symmetry.c get_translation, is_identity = 0).
[[nodiscard]] std::vector<Vector3d>
translations_for_rotation(Cell const &cell, OverlapChecker const &checker,
                          Matrix3i const &rot, int min_index, double symprec) {
  Vector3d const origin = rot.cast<double>() * cell.position(min_index);
  std::vector<Vector3d> result;
  int const n = static_cast<int>(cell.size());
  for (int i = 0; i < n; ++i) {
    if (cell.type(i) != cell.type(min_index)) {
      continue;
    }
    Vector3d const trans = cell.position(i) - origin;
    if (checker.check_total_overlap(trans, rot, symprec))
      result.push_back(math::mod1(trans));
  }
  return result;
}

} // namespace

Result<PointSymmetry> lattice_symmetry(Cell const &cell, double symprec,
                                       AngleTolerance angle_tolerance) {
  auto const min_lattice = reduce::delaunay_reduce(cell.lattice(), symprec);
  if (!min_lattice) {
    return leaf::new_error(e_symmetry_operation_search_failed{});
  }
  Matrix3d const metric_orig = math::metric_tensor(*min_lattice);
  AngleTolerance angle = angle_tolerance;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    auto found =
        collect_metric_symmetries(*min_lattice, metric_orig, symprec, angle);
    if (found) {
      return transform_pointsymmetry(*found, cell.lattice(), *min_lattice);
    }
    if (angle) {
      *angle *= kAngleReduceRate; // too many: tighten and retry
    }
    // std::nullopt (auto): identical retry, will exhaust attempts and fail.
  }
  return leaf::new_error(e_symmetry_operation_search_failed{});
}

Result<SymmetryOperations> find_symmetry(Cell const &cell, double symprec,
                                         AngleTolerance angle_tolerance) {
  BOOST_LEAF_AUTO(lat_sym, lattice_symmetry(cell, symprec, angle_tolerance));
  if (lat_sym.empty())
    return leaf::new_error(e_symmetry_operation_search_failed{});

  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return leaf::new_error(e_empty_cell{});

  OverlapChecker const checker(cell);
  SymmetryOperations ops;
  for (Matrix3i const &rot : lat_sym)
    for (Vector3d const &t :
         translations_for_rotation(cell, checker, rot, *min_index, symprec))
      ops.push_back({rot, t});

  return ops;
}

std::vector<Vector3d> pure_translations(Cell const &cell, double symprec) {
  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return {};
  OverlapChecker const checker(cell);
  return translations_for_rotation(cell, checker, Matrix3i::Identity(),
                                   *min_index, symprec);
}

} // namespace spglib::symmetry
