#include <spglib/core/overlap.hpp>

#include <spglib/math/integer_matrix.hpp>

#include <algorithm>
#include <numeric>
#include <ranges>

namespace spglib {

namespace {

[[nodiscard]] Vector3d row(Positions const &p, int i) noexcept {
  return p.row(i).transpose();
}

// Squared Cartesian distance from a fractional position to its nearest lattice
// point. This key is invariant under integer lattice translations, so it is
// stable under the (mod 1) translation part of a symmetry operation.
[[nodiscard]] double lattice_point_distance_sq(Matrix3d const &lattice,
                                               Vector3d const &frac) noexcept {
  return (lattice * frac.unaryExpr([](double x) { return x - math::nint(x); }))
      .squaredNorm();
}

// Permutation sorting atoms by (type, distance-to-nearest-lattice-point), the
// ordering spglib uses to cluster coincident atoms (overlap.c perm_argsort).
[[nodiscard]] std::vector<int>
argsort_by_distance(Matrix3d const &lattice, Positions const &pos,
                    std::vector<int> const &types) {
  int const n = static_cast<int>(types.size());
  std::vector<double> dist(static_cast<std::size_t>(n));
  std::ranges::transform(std::views::iota(0, n), dist.begin(), [&](int i) {
    return lattice_point_distance_sq(lattice, row(pos, i));
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
                double symprec) noexcept {
  Vector3d diff = (a - b).unaryExpr([](double x) { return x - math::nint(x); });
  return (lattice * diff).norm() <= symprec;
}

OverlapChecker::OverlapChecker(Cell const &cell)
    : lattice_(cell.lattice()), pos_sorted_(cell.size(), 3),
      types_sorted_(static_cast<std::size_t>(cell.size())),
      size_(static_cast<int>(cell.size())) {
  auto const perm =
      argsort_by_distance(cell.lattice(), cell.positions(), cell.types());
  for (int i = 0; i < size_; ++i) {
    pos_sorted_.row(i) =
        cell.positions().row(perm[static_cast<std::size_t>(i)]);
    types_sorted_[static_cast<std::size_t>(i)] =
        cell.types()[static_cast<std::size_t>(
            perm[static_cast<std::size_t>(i)])];
  }
}

bool OverlapChecker::possible_overlap(Vector3d const &trans,
                                      Matrix3i const &rot,
                                      double symprec) const {
  int const search_num = std::min(size_, 3);
  Matrix3d const rotd = rot.cast<double>();
  for (int it = 0; it < search_num; ++it) {
    Vector3d const pr = rotd * row(pos_sorted_, it) + trans;
    int const tr = types_sorted_[static_cast<std::size_t>(it)];
    bool found = false;
    for (int i = 0; i < size_; ++i) {
      if (is_overlap_same_type(pr, row(pos_sorted_, i), tr,
                               types_sorted_[static_cast<std::size_t>(i)],
                               lattice_, symprec)) {
        found = true;
        break;
      }
    }
    if (!found)
      return false;
  }
  return true;
}

bool OverlapChecker::check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot,
                                         double symprec) const {
  if (!possible_overlap(trans, rot, symprec))
    return false;

  bool const is_identity = rot == Matrix3i::Identity();
  Matrix3d const rotd = rot.cast<double>();
  Positions rotated(size_, 3);
  for (int i = 0; i < size_; ++i) {
    Vector3d const p = is_identity ? row(pos_sorted_, i)
                                   : Vector3d(rotd * row(pos_sorted_, i));
    rotated.row(i) = (p + trans).transpose();
  }

  // Greedy matching: every original atom must overlap a not-yet-matched rotated
  // atom (overlap.c check_total_overlap_for_sorted).
  std::vector<char> found(static_cast<std::size_t>(size_), 0);
  int search_start = 0;
  for (int io = 0; io < size_; ++io) {
    while (search_start < size_ &&
           found[static_cast<std::size_t>(search_start)])
      ++search_start;
    int ir = search_start;
    for (; ir < size_; ++ir) {
      if (found[static_cast<std::size_t>(ir)])
        continue;
      if (is_overlap_same_type(row(pos_sorted_, io), row(rotated, ir),
                               types_sorted_[static_cast<std::size_t>(io)],
                               types_sorted_[static_cast<std::size_t>(ir)],
                               lattice_, symprec)) {
        found[static_cast<std::size_t>(ir)] = 1;
        break;
      }
    }
    if (ir == size_)
      return false;
  }
  return true;
}

} // namespace spglib
