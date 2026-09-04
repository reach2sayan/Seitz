#include "symmetry/search.hpp"

#include "core/matrix_order.hpp"
#include "core/overlap.hpp"
#include "core/validation.hpp"
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include "math/lattice_parameters.hpp"
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/periodicity.hpp>

#include <boost/container/flat_map.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <optional>
#include <ranges>
#include <utility>

// The symmetry-operation search:
//   1. Delaunay-reduce the lattice and enumerate the lattice point group by
//      trying every unimodular integer matrix built from 26 candidate axes that
//      preserves the metric.
//   2. Transform those rotations into the input cell's basis.
//   3. For each rotation, find the translations that map the cell onto itself
//      via the OverlapChecker.
namespace cppcrystal::symmetry {

namespace {

constexpr int kNumAttempt = 100;
constexpr double kAngleReduceRate = 0.95;
constexpr double kSinDtheta2Cutoff = 1e-12;

// 26 candidate lattice vectors.
constexpr int kRelativeAxes[26][3] = {
    {1, 0, 0},   {0, 1, 0},   {0, 0, 1},  {-1, 0, 0},  {0, -1, 0},
    {0, 0, -1},  {0, 1, 1},   {1, 0, 1},  {1, 1, 0},   {0, -1, -1},
    {-1, 0, -1}, {-1, -1, 0}, {0, 1, -1}, {-1, 0, 1},  {1, -1, 0},
    {0, -1, 1},  {1, 0, -1},  {-1, 1, 0}, {1, 1, 1},   {-1, -1, -1},
    {-1, 1, 1},  {1, -1, 1},  {1, 1, -1}, {1, -1, -1}, {-1, 1, -1},
    {-1, -1, 1},
};

[[nodiscard]] Vector3i axis(int const (&a)[3]) { return {a[0], a[1], a[2]}; }

// The candidate axis triples worth testing at all: those whose matrix is
// unimodular. That is a property of kRelativeAxes alone, so it is decided once
// at compile time rather than 26^3 times per pass -- and the search re-enters
// this loop up to kNumAttempt times per determination attempt. Index triples
// rather than the matrices themselves: 3 bytes an entry instead of 9, and the
// unpacking below is the same comma-initializer the runtime loop already used.
//
// Order is the cartesian-product order of the loop this replaces, and the two
// filters (determinant here, couples_aperiodic at runtime) commute, so exactly
// the same candidates are visited in exactly the same sequence.
using AxisTriple = std::array<std::uint8_t, 3>;

// 6960 of the 17576 triples are unimodular; `count` below pins that.
struct UnimodularAxes {
  std::array<AxisTriple, 6960> triples{};
  std::size_t count = 0;
};

inline constexpr UnimodularAxes kUnimodularAxes = [] {
  UnimodularAxes out{};
  constexpr std::size_t n = std::size(kRelativeAxes);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      for (std::size_t k = 0; k < n; ++k) {
        // Row-major proxy whose columns are the three candidate vectors:
        // Eigen fixed matrices are not literal types for arithmetic.
        math::Mat3Rows m{};
        for (std::size_t r = 0; r < 3; ++r) {
          m[r * 3 + 0] = kRelativeAxes[i][r];
          m[r * 3 + 1] = kRelativeAxes[j][r];
          m[r * 3 + 2] = kRelativeAxes[k][r];
        }
        if (int const det = math::determinant(m); det == 1 || det == -1) {
          // Overflows at compile time if the size above is ever wrong.
          out.triples[out.count++] = {static_cast<std::uint8_t>(i),
                                      static_cast<std::uint8_t>(j),
                                      static_cast<std::uint8_t>(k)};
        }
      }
    }
  }
  return out;
}();
static_assert(kUnimodularAxes.count == kUnimodularAxes.triples.size());

// Every pairwise dot product among the 26 candidate axes mapped into Cartesian
// space. The metric of a candidate basis (a_i, a_j, a_k) IS the 3x3 of those
// dot products -- metric(p, q) = (L.a_p).(L.a_q) -- so the whole scan shares
// 676 values computed once, instead of each of its 6960 candidates doing two
// 3x3 matrix products to rediscover six of them.
struct AxisMetrics {
  static constexpr std::size_t kAxes = std::size(kRelativeAxes);
  std::array<double, kAxes * kAxes> dot{};
  // sqrt of the diagonal: the length of a candidate basis vector depends only
  // on which of the 26 axes it is, not on the triple it appears in, so the
  // scan reads 26 roots instead of taking 3 per candidate (20880 per pass).
  std::array<double, kAxes> length{};

  [[nodiscard]] double operator()(std::size_t a, std::size_t b) const noexcept {
    return dot[a * kAxes + b];
  }
};

[[nodiscard]] AxisMetrics axis_metrics(Matrix3d const &lattice) {
  std::array<Vector3d, AxisMetrics::kAxes> cartesian;
  for (auto const [a, v] : cartesian | std::views::enumerate) {
    v = lattice * axis(kRelativeAxes[a]).cast<double>();
  }
  AxisMetrics out;
  for (std::size_t a = 0; a < AxisMetrics::kAxes; ++a) {
    for (std::size_t b = 0; b < AxisMetrics::kAxes; ++b) {
      out.dot[a * AxisMetrics::kAxes + b] = cartesian[a].dot(cartesian[b]);
    }
    out.length[a] = std::sqrt(out.dot[a * AxisMetrics::kAxes + a]);
  }
  return out;
}

// The `orig` side of is_identity_metric, which is the same for every candidate
// in a pass: three sqrt and (on the angle path) three acos that used to be
// recomputed 17576 times.
struct MetricReference {
  Eigen::Array3d length;            // sqrt of the metric diagonal
  std::array<double, 3> cosine{};   // pairs (0,1), (0,2), (1,2)
  std::array<double, 3> angle_deg{};
};

// The pairs of basis vectors whose angles the metric comparison tests.
constexpr std::array<std::pair<int, int>, 3> kElemSets{{{0, 1}, {0, 2}, {1, 2}}};

// The degree conversion keeps the original `/ pi * 180` spelling rather than
// folding a 180/pi constant, so the result stays bit-identical to what the
// reference produces.
[[nodiscard]] double metric_angle_deg(Matrix3d const &m, int a, int b) {
  return math::metric_angle(m, a, b) / std::numbers::pi * 180.0;
}

[[nodiscard]] MetricReference reference_of(Matrix3d const &metric) {
  MetricReference ref{.length = metric.diagonal().array().sqrt()};
  for (auto const [i, pair] : kElemSets | std::views::enumerate) {
    auto const [j, k] = pair;
    ref.cosine[static_cast<std::size_t>(i)] = math::metric_cosine(metric, j, k);
    ref.angle_deg[static_cast<std::size_t>(i)] = metric_angle_deg(metric, j, k);
  }
  return ref;
}

// For a layer group the aperiodic axis must map to itself: no candidate axis
// matrix may couple the aperiodic axis with the periodic plane. True if the
// off-diagonal block for `aperiodic_axis` has any nonzero entry. Candidate
// vectors are the columns of `axes`.
[[nodiscard]] bool couples_aperiodic(Matrix3i const &axes,
                                     int aperiodic_axis) noexcept {
  return std::ranges::any_of(std::views::iota(0, 3), [&](int p) {
    return p != aperiodic_axis &&
           (axes(p, aperiodic_axis) != 0 || axes(aperiodic_axis, p) != 0);
  });
}

// Whether a candidate metric agrees with the reference one. An unset
// angle_tolerance uses the sin-based criterion; a positive angle_tolerance
// compares angles in degrees. math::metric_cosine clamps to [-1, 1] -- see
// there for why.
// The cheap half of the comparison: the basis-vector lengths must agree. Split
// out so the scan can run it off three table lookups and skip assembling the
// rest of the metric for the candidates it rejects -- which is most of them.
[[nodiscard]] bool lengths_agree(MetricReference const &orig,
                                 Eigen::Array3d const &length_rot,
                                 double symprec) {
  return (orig.length - length_rot).abs().maxCoeff() <= symprec;
}

// The remaining half: the inter-axis angles must agree too. `length_rot` is
// passed in rather than recomputed, so the three sqrt happen once.
[[nodiscard]] bool angles_agree(Matrix3d const &rotated,
                                Eigen::Array3d const &length_rot,
                                MetricReference const &orig, double symprec,
                                AngleTolerance angle_tolerance) {
  for (auto const [i, pair] : kElemSets | std::views::enumerate) {
    auto const [j, k] = pair;
    auto const ui = static_cast<std::size_t>(i);
    if (angle_tolerance) {
      if (std::fabs(orig.angle_deg[ui] - metric_angle_deg(rotated, j, k)) >
          *angle_tolerance) {
        return false;
      }
    } else {
      double const cos1 = orig.cosine[ui];
      double const cos2 = math::metric_cosine(rotated, j, k);
      double const x =
          cos1 * cos2 + std::sqrt(1 - cos1 * cos1) * std::sqrt(1 - cos2 * cos2);
      double const sin_dtheta2 = 1 - x * x;
      double const length_ave2 =
          ((orig.length[j] + length_rot[j]) * (orig.length[k] + length_rot[k])) /
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
template <GroupFamily F>
[[nodiscard]] std::optional<PointSymmetry> collect_metric_symmetries(
    Matrix3d const &min_lattice, Matrix3d const &metric_orig, double symprec,
    AngleTolerance angle_tolerance, std::optional<int> aperiodic_axis) {
  // Layer groups admit at most 24 lattice symmetries (the aperiodic axis halves
  // the 48 of a 3D point group); overflowing the cap asks for a tighter angle.
  constexpr std::size_t cap = F == GroupFamily::layer ? 24 : 48;
  MetricReference const reference = reference_of(metric_orig);
  AxisMetrics const metrics = axis_metrics(min_lattice);
  PointSymmetry found;
  for (AxisTriple const &triple : kUnimodularAxes.triples) {
    std::array<std::size_t, 3> const a{triple[0], triple[1], triple[2]};

    // The lengths are three lookups, and the length test rejects the large
    // majority of candidates -- so neither the off-diagonals nor the integer
    // axis matrix is assembled until a candidate survives it.
    Eigen::Array3d const length_rot{metrics.length[a[0]], metrics.length[a[1]],
                                    metrics.length[a[2]]};
    if (!lengths_agree(reference, length_rot, symprec)) {
      continue;
    }

    Matrix3i axes;
    axes << axis(kRelativeAxes[a[0]]), axis(kRelativeAxes[a[1]]),
        axis(kRelativeAxes[a[2]]);
    if constexpr (F == GroupFamily::layer) {
      if (aperiodic_axis && couples_aperiodic(axes, *aperiodic_axis)) {
        continue;
      }
    }

    Matrix3d metric;
    for (int p = 0; p < 3; ++p) {
      for (int q = 0; q < 3; ++q) {
        metric(p, q) = metrics(a[static_cast<std::size_t>(p)],
                               a[static_cast<std::size_t>(q)]);
      }
    }
    if (angles_agree(metric, length_rot, reference, symprec,
                     angle_tolerance)) {
      if (found.size() >= cap) {
        return std::nullopt;
      }
      found.push_back(axes);
    }
  }
  return found;
}

// Bring rotations from the Delaunay basis (original_lat) into the cell basis
// (new_lat) by the similarity transform b^-1 . R . b with
// b = original_lat^-1 . new_lat.
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
// empty cell. Among equally rare types the smallest index wins.
[[nodiscard]] std::optional<int> index_with_least_atoms(Cell const &cell) {
  boost::container::flat_map<int, int> frequency;
  for (int const type : cell.types()) {
    ++frequency[type];
  }
  if (frequency.empty()) {
    return std::nullopt;
  }
  int const rarest = std::ranges::min(frequency | std::views::values);
  auto const first = std::ranges::find_if(
      cell.types(), [&](int type) { return frequency.at(type) == rarest; });
  return static_cast<int>(first - cell.types().begin());
}

// Translations t such that x -> rot . x + t maps the cell onto itself.
[[nodiscard]] std::vector<Vector3d>
translations_for_rotation(Cell const &cell, OverlapChecker const &checker,
                          Matrix3i const &rot, int min_index) {
  Vector3d const origin = rot.cast<double>() * cell.position(min_index);
  std::vector<Vector3d> result;
  int const n = static_cast<int>(cell.size());
  for (int i = 0; i < n; ++i) {
    if (cell.type(i) != cell.type(min_index)) {
      continue;
    }
    Vector3d const trans = cell.position(i) - origin;
    if (checker.check_total_overlap(trans, rot)) {
      // Layer translations live in the periodic plane; the aperiodic component
      // is kept raw rather than folded into [0, 1).
      result.push_back(wrap(trans, cell.periodicity()));
    }
  }
  return result;
}

} // namespace

template <GroupFamily F>
Result<PointSymmetry> SymmetrySearch<F>::lattice_symmetry() const {
  Cell const &cell = cell_;
  double const symprec = tol_.symprec;
  std::optional<int> const layer_axis = aperiodic_axis(cell.periodicity());
  // A layer cell reduces only the two periodic lattice vectors, leaving the
  // aperiodic axis fixed; a 3D cell uses the full Delaunay reduction.
  auto const min_lattice = [&] {
    if constexpr (F == GroupFamily::layer) {
      return cell.lattice().delaunay_in_plane(layer_axis.value_or(2), symprec);
    } else {
      return cell.lattice().delaunay(symprec);
    }
  }();
  if (!min_lattice) {
    return leaf::new_error(e_symmetry_operation_search_failed{});
  }
  Matrix3d const metric_orig = min_lattice->metric();
  AngleTolerance angle = tol_.angle_tolerance;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    auto found = collect_metric_symmetries<F>(
        min_lattice->matrix(), metric_orig, symprec, angle, layer_axis);
    if (found) {
      return transform_pointsymmetry(*found, cell.lattice().matrix(),
                                     min_lattice->matrix());
    }
    if (angle) {
      *angle *= kAngleReduceRate; // too many: tighten and retry
    }
    // std::nullopt (auto): identical retry, will exhaust attempts and fail.
  }
  return leaf::new_error(e_symmetry_operation_search_failed{});
}

template <GroupFamily F>
Result<Operations> SymmetrySearch<F>::operations() const {
  Cell const &cell = cell_;
  if (auto valid = validate_cell(cell); !valid) {
    return valid.error();
  }
  BOOST_LEAF_AUTO(lat_sym, lattice_symmetry());
  if (lat_sym.empty())
    return leaf::new_error(e_symmetry_operation_search_failed{});

  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return leaf::new_error(e_empty_cell{});

  OverlapChecker const checker(cell, tol_.symprec);
  std::vector<SymmetryOperation> ops;
  for (Matrix3i const &rot : lat_sym)
    for (Vector3d const &t :
         translations_for_rotation(cell, checker, rot, *min_index))
      ops.push_back({rot, t});

  return Operations{std::move(ops)};
}

template <GroupFamily F>
Operations SymmetrySearch<F>::reduce(Operations const &operations,
                                     Tolerance const &tol) const {
  Cell const &cell = cell_;
  auto const lat_sym = SymmetrySearch{cell, tol}.lattice_symmetry();
  if (!lat_sym) {
    return {};
  }
  OverlapChecker const checker(cell, tol.symprec);
  RotationSet const lattice_rotations = rotation_set(*lat_sym);
  auto survives = [&](SymmetryOperation const &op) {
    return lattice_rotations.contains(op.rotation) &&
           checker.check_total_overlap(op.translation, op.rotation);
  };
  return Operations{std::from_range, operations | std::views::filter(survives)};
}

template <GroupFamily F>
std::vector<Vector3d> SymmetrySearch<F>::pure_translations() const {
  Cell const &cell = cell_;
  std::optional<int> const min_index = index_with_least_atoms(cell);
  if (!min_index)
    return {};
  OverlapChecker const checker(cell, tol_.symprec);
  return translations_for_rotation(cell, checker, Matrix3i::Identity(),
                                   *min_index);
}

template class SymmetrySearch<GroupFamily::space>;
template class SymmetrySearch<GroupFamily::layer>;

} // namespace cppcrystal::symmetry
