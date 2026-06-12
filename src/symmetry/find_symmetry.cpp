#include <spglib/symmetry/find_symmetry.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/core/validation.hpp>
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

// For a layer group the aperiodic axis must map to itself: no candidate axis
// matrix may couple the aperiodic axis with the periodic plane. True if the
// off-diagonal block for `aperiodic_axis` has any nonzero entry (symmetry.c's
// per-axis switch). Candidate vectors are the columns of `axes`, so axes(r,c)
// matches the reference axes[r][c].
[[nodiscard]] bool couples_aperiodic(Matrix3i const &axes,
                                     int aperiodic_axis) noexcept {
  return std::ranges::any_of(std::views::iota(0, 3), [&](int p) {
    return p != aperiodic_axis &&
           (axes(p, aperiodic_axis) != 0 || axes(aperiodic_axis, p) != 0);
  });
}

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
                          AngleTolerance angle_tolerance,
                          std::optional<int> aperiodic_axis) {
  // Layer groups admit at most 24 lattice symmetries (the aperiodic axis halves
  // the 48 of a 3D point group); overflowing the cap asks for a tighter angle.
  std::size_t const cap = aperiodic_axis ? 24 : 48;
  PointSymmetry found;
  for (auto const &ai : kRelativeAxes) {
    for (auto const &aj : kRelativeAxes) {
      for (auto const &ak : kRelativeAxes) {
        Matrix3i axes;
        axes << axis(ai), axis(aj), axis(ak);
        if (aperiodic_axis && couples_aperiodic(axes, *aperiodic_axis)) {
          continue;
        }
        if (int const det = axes.determinant(); det != 1 && det != -1) {
          continue;
        }
        Matrix3d const lattice = min_lattice * axes.cast<double>();
        Matrix3d const metric = math::metric_tensor(lattice);
        if (is_identity_metric(metric, metric_orig, symprec, angle_tolerance)) {
          if (found.size() >= cap) {
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
                          Matrix3i const &rot, int min_index, double symprec,
                          std::optional<int> aperiodic_axis) {
  Vector3d const origin = rot.cast<double>() * cell.position(min_index);
  std::vector<Vector3d> result;
  int const n = static_cast<int>(cell.size());
  for (int i = 0; i < n; ++i) {
    if (cell.type(i) != cell.type(min_index)) {
      continue;
    }
    Vector3d const trans = cell.position(i) - origin;
    if (checker.check_total_overlap(trans, rot, symprec)) {
      // Layer translations live in the periodic plane; the aperiodic component
      // is kept raw rather than folded into [0, 1) (symmetry.c get_layer_translation).
      Vector3d wrapped;
      for (int j = 0; j < 3; ++j) {
        wrapped[j] = (aperiodic_axis && j == *aperiodic_axis)
                         ? trans[j]
                         : math::wrap_to_unit_cell(trans[j]);
      }
      result.push_back(wrapped);
    }
  }
  return result;
}

} // namespace

Result<PointSymmetry> lattice_symmetry(Cell const &cell, double symprec,
                                       AngleTolerance angle_tolerance) {
  std::optional<int> const aperiodic_axis = cell.aperiodic_axis();
  // A layer cell reduces only the two periodic lattice vectors, leaving the
  // aperiodic axis fixed (del_layer_delaunay_reduce); a 3D cell uses the full
  // Delaunay reduction.
  auto const min_lattice =
      aperiodic_axis
          ? reduce::delaunay_reduce(cell.lattice(), *aperiodic_axis, symprec)
          : reduce::delaunay_reduce(cell.lattice(), symprec);
  if (!min_lattice) {
    return leaf::new_error(e_symmetry_operation_search_failed{});
  }
  Matrix3d const metric_orig = math::metric_tensor(*min_lattice);
  AngleTolerance angle = angle_tolerance;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    auto found = collect_metric_symmetries(*min_lattice, metric_orig, symprec,
                                           angle, aperiodic_axis);
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
  if (auto valid = validate_cell(cell); !valid) {
    return valid.error();
  }
  BOOST_LEAF_AUTO(lat_sym, lattice_symmetry(cell, symprec, angle_tolerance));
  if (lat_sym.empty())
    return leaf::new_error(e_symmetry_operation_search_failed{});

  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return leaf::new_error(e_empty_cell{});

  OverlapChecker const checker(cell);
  std::optional<int> const aperiodic_axis = cell.aperiodic_axis();
  SymmetryOperations ops;
  for (Matrix3i const &rot : lat_sym)
    for (Vector3d const &t : translations_for_rotation(
             cell, checker, rot, *min_index, symprec, aperiodic_axis))
      ops.push_back({rot, t});

  return ops;
}

SymmetryOperations reduce_symmetry(Cell const &cell,
                                   SymmetryOperations const &operations,
                                   double symprec,
                                   AngleTolerance angle_tolerance) {
  auto const lat_sym = lattice_symmetry(cell, symprec, angle_tolerance);
  if (!lat_sym) {
    return {};
  }
  OverlapChecker const checker(cell);
  auto survives = [&](SymmetryOperation const &op) {
    return std::ranges::any_of(
               *lat_sym, [&](Matrix3i const &r) { return r == op.rotation; }) &&
           checker.check_total_overlap(op.translation, op.rotation, symprec);
  };
  auto kept = operations | std::views::filter(survives);
  return {kept.begin(), kept.end()};
}

std::vector<Vector3d> pure_translations(Cell const &cell, double symprec) {
  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return {};
  OverlapChecker const checker(cell);
  return translations_for_rotation(cell, checker, Matrix3i::Identity(),
                                   *min_index, symprec, cell.aperiodic_axis());
}

} // namespace spglib::symmetry
