#include <cppcrystal/core/position_index.hpp>

#include <cppcrystal/math/fractional.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <ranges>

namespace cppcrystal {

namespace {

// Upper bound on buckets per periodic axis: beyond this a bucket is already
// far narrower than any coordinate noise, and the key stays well inside int32.
constexpr double kMaxDivisions = double{1 << 20};

[[nodiscard]] std::int64_t floor_to_int64(double x) noexcept {
  return static_cast<std::int64_t>(std::floor(x));
}

} // namespace

bool coincident(Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
                double symprec, CellPeriodicity const &periodicity) noexcept {
  return (lattice * minimal_image(a - b, periodicity)).norm() <= symprec;
}

BucketGeometry BucketGeometry::of(Matrix3d const &lattice, double symprec,
                                  CellPeriodicity const &periodicity) noexcept {
  Matrix3d const inv = lattice.inverse();
  BucketGeometry g{};
  for (auto const [axis, kind] : periodicity | std::views::enumerate) {
    auto const i = static_cast<Index>(axis);
    // Twice the largest fractional displacement a symprec-sized Cartesian
    // step can cause along this axis.
    double const bound = 2.0 * symprec * inv.row(i).norm();
    if (kind == AxisKind::periodic) {
      double const wanted = bound > 0.0 ? 1.0 / bound : kMaxDivisions;
      auto const divisions =
          static_cast<int>(std::floor(std::clamp(wanted, 1.0, kMaxDivisions)));
      g.divisions[i] = divisions;
      g.width[i] = 1.0 / divisions;
    } else {
      g.divisions[i] = 0;
      g.width[i] =
          bound > 0.0 ? bound : std::numeric_limits<double>::infinity();
    }
  }
  return g;
}

BucketGeometry::Key
BucketGeometry::bucket_of(Vector3d const &frac) const noexcept {
  Key key{};
  for (Index i = 0; i < 3; ++i) {
    if (divisions[i] > 0) {
      auto const raw =
          floor_to_int64(math::wrap_to_unit_cell(frac[i]) * divisions[i]);
      // wrap_to_unit_cell lands in [0, 1), but f * div can round up to div.
      key[static_cast<std::size_t>(i)] =
          std::min(raw, std::int64_t{divisions[i]} - 1);
    } else {
      key[static_cast<std::size_t>(i)] = floor_to_int64(frac[i] / width[i]);
    }
  }
  return key;
}

PositionIndex::PositionIndex(BucketGeometry geometry, Positions const &positions,
                             Types const &types, Matrix3d const &lattice,
                             double symprec, CellPeriodicity const &periodicity)
    : geometry_(geometry), lattice_(lattice), symprec_(symprec),
      periodicity_(periodicity), positions_(&positions), types_(&types) {
  entries_.reserve(types.size());
  for (Index i = 0; i < positions.rows(); ++i) {
    entries_.push_back(
        {geometry_.bucket_of(positions.row(i).transpose()), static_cast<int>(i)});
  }
  std::ranges::sort(entries_);
}

PositionIndex::PositionIndex(Cell const &cell, double symprec)
    : PositionIndex(BucketGeometry::of(cell.lattice().matrix(), symprec,
                                       cell.periodicity()),
                    cell.positions(), cell.types(), cell.lattice().matrix(), symprec,
                    cell.periodicity()) {}

PositionIndex::Buckets
PositionIndex::bucket_ranges(Vector3d const &point) const {
  auto const centre = geometry_.bucket_of(point);
  auto const neighbour = [&](std::size_t axis, int offset) {
    auto const b = centre[axis] + offset;
    std::int64_t const div = geometry_.divisions[static_cast<Index>(axis)];
    return div > 0 ? (b + div) % div : b;
  };

  // With one or two divisions on an axis the neighbourhood folds onto itself;
  // sort + unique visits every bucket exactly once.
  boost::container::static_vector<BucketGeometry::Key, 27> keys;
  auto const offsets = std::views::iota(-1, 2);
  for (auto const [d0, d1, d2] :
       std::views::cartesian_product(offsets, offsets, offsets)) {
    keys.push_back({neighbour(0, d0), neighbour(1, d1), neighbour(2, d2)});
  }
  std::ranges::sort(keys);
  auto const [dup, end] = std::ranges::unique(keys);
  keys.erase(dup, end);

  Buckets out;
  for (auto const &key : keys) {
    auto const run = std::ranges::equal_range(entries_, key, {}, &Entry::bucket);
    if (!run.empty()) {
      out.push_back(std::span<Entry const>(run.begin(), run.end()));
    }
  }
  return out;
}

bool PositionIndex::coincides(Vector3d const &point, int atom) const noexcept {
  return coincident(point, positions_->row(atom).transpose(), lattice_,
                    symprec_, periodicity_);
}

} // namespace cppcrystal
