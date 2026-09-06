#include "core/position_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <numeric>
#include <ranges>

namespace seitz {

namespace {

// Safety factor on the derived image reach. The bound below is exact, but it
// is cheap to be generous -- an over-estimate only scans a bucket that turns
// out to be empty, while an under-estimate would drop a real match.
constexpr double kImageMargin = 2.0;

// The largest fractional offset a Cartesian offset of `half_width` can be:
// ||lattice^-1||_inf * half_width
[[nodiscard]] double image_reach_of(Matrix3d const &lattice_inv,
                                    double half_width) noexcept {
  double const reach = kImageMargin *
                       lattice_inv.cwiseAbs().rowwise().sum().maxCoeff() *
                       half_width;
  return std::isfinite(reach) ? reach : kImageMargin;
}

// Buckets per periodic axis: about one atom per bucket, but never so many
// that a query's reach spans more than one of them anyway -- a generation
// distance check, whose "symprec" is a bond length, gets a single bucket and
// therefore the plain pairwise scan, with no bucketing on top.
[[nodiscard]] int divisions_for(Index n, double reach) noexcept {
  double const by_atoms = std::cbrt(static_cast<double>(n));
  double const by_reach = reach > 0.0 ? 0.5 / reach : 64.0;
  return std::clamp(static_cast<int>(std::min(by_atoms, by_reach)), 1, 64);
}

[[nodiscard]] int wrap_index(long i, int n) noexcept {
  long const m = i % n;
  return static_cast<int>(m < 0 ? m + n : m);
}

} // namespace

PositionIndex::PositionIndex(Positions const &positions, Types const &types,
                             Matrix3d const &lattice, double symprec,
                             CellPeriodicity const &periodicity)
    : lattice_(lattice), symprec_(symprec),
      // The sphere test in `coincident` and the bucket walk here round
      // differently; a few ulps of the coordinate scale keep a match on the
      // sphere's surface from landing just outside the scanned buckets.
      image_reach_(image_reach_of(
          lattice.inverse(),
          symprec + 64 * std::numeric_limits<double>::epsilon() *
                        (symprec + lattice.cwiseAbs().maxCoeff()))),
      periodicity_(periodicity), positions_(&positions), types_(&types) {
  auto const n = static_cast<int>(positions.rows());
  int const divisions = divisions_for(n, image_reach_);
  std::ranges::transform(periodicity, divisions_.begin(), [&](AxisKind kind) {
    return kind == AxisKind::periodic ? divisions : 1;
  });
  everything_ = {std::from_range, std::views::iota(0, n)};
  if (divisions_ == std::array{1, 1, 1}) {
    return; // one bucket: every query is `everything_`
  }
  auto const buckets = static_cast<std::size_t>(
      std::ranges::fold_left(divisions_, 1, std::multiplies{}));

  // Counting sort by bucket: a histogram, its prefix sums as bucket starts,
  // and a stable scatter that keeps each bucket in index order.
  std::vector<int> const bucket_of{
      std::from_range, std::views::iota(0, n) | std::views::transform([&](int i) {
                         return bucket(wrap(positions.row(i).transpose(),
                                            periodicity));
                       })};
  starts_.assign(buckets + 1, 0);
  for (int const b : bucket_of) {
    ++starts_[static_cast<std::size_t>(b) + 1];
  }
  std::inclusive_scan(starts_.begin(), starts_.end(), starts_.begin());
  std::vector<int> cursor(starts_.begin(), std::prev(starts_.end()));
  atoms_.resize(static_cast<std::size_t>(n));
  for (auto const [i, b] : std::views::enumerate(bucket_of)) {
    atoms_[static_cast<std::size_t>(cursor[static_cast<std::size_t>(b)]++)] =
        static_cast<int>(i);
  }
}

int PositionIndex::bucket(Vector3d const &folded) const noexcept {
  return (bucket_along(folded[0], 0) * divisions_[1] +
          bucket_along(folded[1], 1)) *
             divisions_[2] +
         bucket_along(folded[2], 2);
}

PositionIndex::PositionIndex(Cell const &cell, double symprec)
    : PositionIndex(cell.positions(), cell.types(), cell.lattice().matrix(),
                    symprec, cell.periodicity()) {}

int PositionIndex::bucket_along(double x, std::size_t axis) const noexcept {
  int const n = divisions_[axis];
  return n == 1 ? 0 : wrap_index(static_cast<long>(std::floor(x * n)), n);
}

std::span<int const> PositionIndex::candidates(Vector3d const &point,
                                               Scratch &out) const {
  if (divisions_ == std::array{1, 1, 1}) {
    return everything_;
  }
  out.clear();
  Vector3d const folded = wrap(point, periodicity_);

  // Per axis, the cyclic run of buckets covering [x - reach, x + reach]: its
  // first bucket and how many. A run that would lap the axis is the whole
  // axis. Bucket indices are taken before wrapping so the count is exact.
  std::array<int, 3> first{};
  std::array<int, 3> count{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    int const n = divisions_[axis];
    if (n == 1) {
      count[axis] = 1;
      continue;
    }
    double const x = folded[static_cast<Index>(axis)];
    auto const lo = static_cast<long>(std::floor((x - image_reach_) * n));
    auto const hi = static_cast<long>(std::floor((x + image_reach_) * n));
    if (hi - lo + 1 >= n) {
      count[axis] = n;
    } else {
      first[axis] = wrap_index(lo, n);
      count[axis] = static_cast<int>(hi - lo + 1);
    }
  }

  if (count == divisions_) {
    return everything_; // already ascending, nothing to gather
  }
  auto const run = [&](std::size_t axis) {
    return std::views::iota(0, count[axis]) |
           std::views::transform([=, this](int i) {
             return wrap_index(first[axis] + i, divisions_[axis]);
           });
  };
  for (auto const [b0, b1, b2] :
       std::views::cartesian_product(run(0), run(1), run(2))) {
    auto const b =
        static_cast<std::size_t>((b0 * divisions_[1] + b1) * divisions_[2] + b2);
    out.insert(out.end(), atoms_.begin() + starts_[b],
               atoms_.begin() + starts_[b + 1]);
  }

  // Each bucket is ascending and holds each atom once; several buckets need
  // one sort to restore the ascending contract.
  if (count[0] * count[1] * count[2] > 1) {
    std::ranges::sort(out);
  }
  return {out.data(), out.size()};
}

std::vector<int> PositionIndex::candidates(Vector3d const &point) const {
  Scratch scratch;
  auto const found = candidates(point, scratch);
  return {found.begin(), found.end()};
}

} // namespace seitz
