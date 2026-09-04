#include "core/position_index.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>

namespace cppcrystal {

namespace bgi = boost::geometry::index;

namespace {

// Safety factor on the derived image reach. The bound below is exact, but it
// is cheap to be generous -- an over-estimate only queries a box that turns
// out to be empty, while an under-estimate would drop a real match.
constexpr double kImageMargin = 2.0;

// The largest fractional offset a Cartesian offset of `half_width` can be:
// ||lattice^-1||_inf * half_width, the induced infinity norm being the maximum
// absolute row sum. A singular (or near-singular) lattice yields a non-finite
// bound; fall back to a reach that admits every image, which is what the
// unpruned search did.
[[nodiscard]] double image_reach_of(Matrix3d const &lattice_inv,
                                    double half_width) noexcept {
  double const reach = kImageMargin *
                       lattice_inv.cwiseAbs().rowwise().sum().maxCoeff() *
                       half_width;
  return std::isfinite(reach) ? reach : 2.0;
}

} // namespace

bool coincident(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
                double symprec, CellPeriodicity const &periodicity) noexcept {
  return (lattice * minimal_image(a - b, periodicity)).norm() <= symprec;
}

PositionIndex::PositionIndex(Positions const &positions, Types const &types,
                             Matrix3d const &lattice, double symprec,
                             CellPeriodicity const &periodicity)
    : lattice_(lattice), lattice_inv_(lattice.inverse()), symprec_(symprec),
      // The sphere test in `coincident` and the box test here round
      // differently; a few ulps of the coordinate scale keep a match on the
      // sphere's surface from landing just outside the box.
      half_width_(symprec + 64 * std::numeric_limits<double>::epsilon() *
                                (symprec + lattice.cwiseAbs().maxCoeff())),
      image_reach_(image_reach_of(lattice_inv_, half_width_)),
      periodicity_(periodicity), positions_(&positions), types_(&types) {
  std::vector<Tree::value_type> values;
  values.reserve(types.size());
  for (Index i = 0; i < positions.rows(); ++i) {
    Vector3d const c = cartesian(positions.row(i).transpose());
    values.emplace_back(Point(c[0], c[1], c[2]), static_cast<int>(i));
  }
  tree_ = Tree(values);
}

PositionIndex::PositionIndex(Cell const &cell, double symprec)
    : PositionIndex(cell.positions(), cell.types(), cell.lattice().matrix(),
                    symprec, cell.periodicity()) {}

Vector3d PositionIndex::cartesian(Vector3d const &frac) const noexcept {
  return lattice_ * wrap(frac, periodicity_);
}

std::span<int const> PositionIndex::candidates(Vector3d const &point,
                                               Scratch &out) const {
  out.clear();

  Vector3d const folded = wrap(point, periodicity_);
  Vector3d const base = lattice_ * folded;
  double const h = half_width_;

  // The images worth querying, per axis. Both `folded` and the indexed atoms
  // live in [0, 1), so the +1 image can only hold a match when the folded
  // coordinate is within image_reach_ of 0, and the -1 image only when it is
  // within image_reach_ of 1. An aperiodic axis has no images at all.
  std::array<std::array<int, 3>, 3> shifts{};
  std::array<int, 3> counts{};
  for (std::size_t axis = 0; axis < 3; ++axis) {
    int n = 0;
    shifts[axis][static_cast<std::size_t>(n++)] = 0;
    if (periodicity_[axis] == AxisKind::periodic) {
      if (folded[static_cast<Index>(axis)] < image_reach_) {
        shifts[axis][static_cast<std::size_t>(n++)] = 1;
      }
      if (folded[static_cast<Index>(axis)] > 1.0 - image_reach_) {
        shifts[axis][static_cast<std::size_t>(n++)] = -1;
      }
    }
    counts[axis] = n;
  }
  bool const single = counts[0] == 1 && counts[1] == 1 && counts[2] == 1;

  for (int i0 = 0; i0 < counts[0]; ++i0) {
    for (int i1 = 0; i1 < counts[1]; ++i1) {
      for (int i2 = 0; i2 < counts[2]; ++i2) {
        int const s0 = shifts[0][static_cast<std::size_t>(i0)];
        int const s1 = shifts[1][static_cast<std::size_t>(i1)];
        int const s2 = shifts[2][static_cast<std::size_t>(i2)];
        // base + lattice * (s0, s1, s2), written as column combinations
        // because the shifts are 0 or +-1.
        Vector3d c = base;
        if (s0 != 0) {
          c += static_cast<double>(s0) * lattice_.col(0);
        }
        if (s1 != 0) {
          c += static_cast<double>(s1) * lattice_.col(1);
        }
        if (s2 != 0) {
          c += static_cast<double>(s2) * lattice_.col(2);
        }
        Box const box(Point(c[0] - h, c[1] - h, c[2] - h),
                      Point(c[0] + h, c[1] + h, c[2] + h));
        for (auto it = tree_.qbegin(bgi::intersects(box)); it != tree_.qend();
             ++it) {
          out.push_back(it->second);
        }
      }
    }
  }

  // The tree yields each atom once per box, so duplicates only arise when a
  // second image was queried -- and the ascending order is contract.
  std::ranges::sort(out);
  if (!single) {
    auto const [dup, end] = std::ranges::unique(out);
    out.erase(dup, end);
  }
  return {out.data(), out.size()};
}

std::vector<int> PositionIndex::candidates(Vector3d const &point) const {
  Scratch scratch;
  auto const found = candidates(point, scratch);
  return {found.begin(), found.end()};
}

bool PositionIndex::coincides(Vector3d const &point, int atom) const noexcept {
  return coincident(point, positions_->row(atom).transpose(), lattice_,
                    symprec_, periodicity_);
}

} // namespace cppcrystal
