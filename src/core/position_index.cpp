#include "core/position_index.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <ranges>
#include <span>

namespace cppcrystal {

namespace bgi = boost::geometry::index;

bool coincident(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
                double symprec, CellPeriodicity const &periodicity) noexcept {
  return (lattice * minimal_image(a - b, periodicity)).norm() <= symprec;
}

PositionIndex::PositionIndex(Positions const &positions, Types const &types,
                             Matrix3d const &lattice, double symprec,
                             CellPeriodicity const &periodicity)
    : lattice_(lattice), symprec_(symprec),
      // The sphere test in `coincident` and the box test here round
      // differently; a few ulps of the coordinate scale keep a match on the
      // sphere's surface from landing just outside the box.
      half_width_(symprec + 64 * std::numeric_limits<double>::epsilon() *
                                (symprec + lattice.cwiseAbs().maxCoeff())),
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

std::vector<int> PositionIndex::candidates(Vector3d const &point) const {
  static constexpr std::array<int, 3> kImages{-1, 0, 1};
  static constexpr std::array<int, 1> kNone{0};
  auto const shifts = [&](std::size_t axis) {
    return periodicity_[axis] == AxisKind::periodic
               ? std::span<int const>(kImages)
               : std::span<int const>(kNone);
  };

  Vector3d const folded = wrap(point, periodicity_);
  double const h = half_width_;
  std::vector<int> out;
  for (auto const [s0, s1, s2] :
       std::views::cartesian_product(shifts(0), shifts(1), shifts(2))) {
    Vector3d const c =
        lattice_ * (folded + Vector3i(s0, s1, s2).cast<double>());
    Box const box(Point(c[0] - h, c[1] - h, c[2] - h),
                  Point(c[0] + h, c[1] + h, c[2] + h));
    for (auto it = tree_.qbegin(bgi::intersects(box)); it != tree_.qend();
         ++it) {
      out.push_back(it->second);
    }
  }
  // A box wider than the cell along an axis sees the same atom from two
  // images.
  std::ranges::sort(out);
  auto const [dup, end] = std::ranges::unique(out);
  out.erase(dup, end);
  return out;
}

bool PositionIndex::coincides(Vector3d const &point, int atom) const noexcept {
  return coincident(point, positions_->row(atom).transpose(), lattice_,
                    symprec_, periodicity_);
}

} // namespace cppcrystal
