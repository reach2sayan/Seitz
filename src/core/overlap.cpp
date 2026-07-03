#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <algorithm>
#include <numeric>
#include <ranges>

namespace cppcrystal {

namespace {

[[nodiscard]] Vector3d row(Positions const &p, int i) noexcept {
  return p.row(i).transpose();
}

// Minimal-image fractional offset of `diff`. For a layer cell the aperiodic
// axis is not periodic, so its component is left as the raw difference; only
// the two periodic components are folded.
[[nodiscard]] Vector3d minimal_image(Vector3d const &diff,
                                     std::optional<int> aperiodic_axis) noexcept {
  Vector3d off = math::nearest_offset(diff);
  if (aperiodic_axis) {
    off[*aperiodic_axis] = diff[*aperiodic_axis];
  }
  return off;
}

// Squared Cartesian distance from a fractional position to its nearest lattice
// point. This key is invariant under integer lattice translations, so it is
// stable under the (mod 1) translation part of a symmetry operation. For a
// layer cell only the periodic directions span lattice points.
[[nodiscard]] double
lattice_point_distance_sq(Matrix3d const &lattice, Vector3d const &frac,
                          std::optional<int> aperiodic_axis) noexcept {
  return (lattice * minimal_image(frac, aperiodic_axis)).squaredNorm();
}

// Permutation sorting atoms by (type, distance-to-nearest-lattice-point), the
// ordering used to cluster coincident atoms.
[[nodiscard]] std::vector<int>
argsort_by_distance(Matrix3d const &lattice, Positions const &pos,
                    std::vector<int> const &types,
                    std::optional<int> aperiodic_axis) {
  int const n = static_cast<int>(types.size());
  std::vector<double> dist(static_cast<std::size_t>(n));
  std::ranges::transform(std::views::iota(0, n), dist.begin(), [&](int i) {
    return lattice_point_distance_sq(lattice, row(pos, i), aperiodic_axis);
  });

  std::vector<int> perm(static_cast<std::size_t>(n));
  std::ranges::iota(perm, 0);
  std::ranges::sort(perm, [&](int a, int b) {
    const auto ua = static_cast<std::size_t>(a);
    const auto ub = static_cast<std::size_t>(b);
    return std::tie(types[ua], dist[ua]) < std::tie(types[ub], dist[ub]);
  });

  return perm;
}

} // namespace

bool is_overlap(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
                double symprec, std::optional<int> aperiodic_axis) noexcept {
  Vector3d const diff = minimal_image(a - b, aperiodic_axis);
  return (lattice * diff).norm() <= symprec;
}

OverlapChecker::OverlapChecker(Cell const &cell)
    : lattice_(cell.lattice()), pos_sorted_(cell.size(), 3),
      types_sorted_(static_cast<std::size_t>(cell.size())),
      size_(static_cast<int>(cell.size())),
      aperiodic_axis_(cell.aperiodic_axis()) {
  auto const perm = argsort_by_distance(cell.lattice(), cell.positions(),
                                        cell.types(), aperiodic_axis_);
  for (auto const [i, src] : perm | std::views::enumerate) {
    pos_sorted_.row(i) = cell.positions().row(src);
    types_sorted_[static_cast<std::size_t>(i)] =
        cell.types()[static_cast<std::size_t>(src)];
  }
}

bool OverlapChecker::possible_overlap(Vector3d const &trans,
                                      Matrix3i const &rot,
                                      double symprec) const {
  int const search_num = std::min(size_, 3);
  Matrix3d const rotd = rot.cast<double>();
  // Every probed atom must coincide with some atom of the same type.
  return std::ranges::all_of(std::views::iota(0, search_num), [&](int it) {
    Vector3d const pr = rotd * row(pos_sorted_, it) + trans;
    int const tr = types_sorted_[static_cast<std::size_t>(it)];
    return std::ranges::any_of(std::views::iota(0, size_), [&](int i) {
      return is_overlap_same_type(pr, row(pos_sorted_, i), tr,
                                  types_sorted_[static_cast<std::size_t>(i)],
                                  lattice_, symprec, aperiodic_axis_);
    });
  });
}

bool OverlapChecker::check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot,
                                         double symprec) const {
  if (!possible_overlap(trans, rot, symprec)) {
    return false;
  }

  bool const is_identity = rot == Matrix3i::Identity();
  Matrix3d const rotd = rot.cast<double>();
  Positions rotated(size_, 3);
  for (int i = 0; i < size_; ++i) {
    Vector3d const p = is_identity ? row(pos_sorted_, i)
                                   : Vector3d(rotd * row(pos_sorted_, i));
    rotated.row(i) = (p + trans).transpose();
  }

  // Greedy matching: every original atom must overlap a not-yet-matched rotated
  // atom.
  std::vector<char> found(static_cast<std::size_t>(size_), 0);
  int search_start = 0;
  for (int io = 0; io < size_; ++io) {
    const auto it =
        std::ranges::find(found.begin() + search_start, found.end(), false);

    search_start = static_cast<int>(std::ranges::distance(found.begin(), it));
    int ir = search_start;
    for (; ir < size_; ++ir) {
      if (found[static_cast<std::size_t>(ir)]) {
        continue;
      }
      if (is_overlap_same_type(row(pos_sorted_, io), row(rotated, ir),
                               types_sorted_[static_cast<std::size_t>(io)],
                               types_sorted_[static_cast<std::size_t>(ir)],
                               lattice_, symprec, aperiodic_axis_)) {
        found[static_cast<std::size_t>(ir)] = 1;
        break;
      }
    }
    if (ir == size_) {
      return false;
    }
  }
  return true;
}

} // namespace cppcrystal
