#include "core/overlap.hpp"

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <ranges>
#include <tuple>
#include <vector>

namespace cppcrystal {

namespace {

// The cell with atoms sorted by (type, squared Cartesian distance to the
// nearest lattice point). The distance is invariant under lattice
// translations, so it is a stable key for clustering coincident atoms; the
// first few atoms in this order are the cheap-rejection probes.
[[nodiscard]] Cell sorted_by_distance(Cell const &cell) {
  struct Keyed {
    int type;
    double distance_sq;
    Index atom;
  };
  std::vector<Keyed> keyed;
  keyed.reserve(static_cast<std::size_t>(cell.size()));
  for (auto const [i, atom] : cell.atoms() | std::views::enumerate) {
    auto const &[position, type] = atom;
    Vector3d const image = minimal_image(position, cell.periodicity());
    keyed.push_back({type, (cell.lattice().matrix() * image).squaredNorm(),
                     static_cast<Index>(i)});
  }
  std::ranges::sort(keyed, {}, [](Keyed const &k) {
    return std::tie(k.type, k.distance_sq);
  });

  Positions positions(cell.size(), 3);
  Types types;
  types.reserve(keyed.size());
  for (auto const [row, k] : keyed | std::views::enumerate) {
    positions.row(row) = cell.positions().row(k.atom);
    types.push_back(k.type);
  }
  return {cell.lattice(), std::move(positions), std::move(types),
          cell.periodicity()};
}

// x -> rot . x + trans applied to every row.
[[nodiscard]] Positions transformed(Positions const &positions,
                                    Matrix3i const &rot,
                                    Vector3d const &trans) {
  Positions out = positions * rot.cast<double>().transpose();
  out.rowwise() += trans.transpose();
  return out;
}

} // namespace

OverlapChecker::OverlapChecker(Cell const &cell, double symprec)
    : sorted_{sorted_by_distance(cell)}, symprec_{symprec},
      index_{sorted_, symprec} {}

bool OverlapChecker::possible_overlap(Positions const &rotated) const {
  Index const probes = std::min<Index>(sorted_.size(), 3);
  return std::ranges::all_of(std::views::iota(Index{0}, probes), [&](Index i) {
    return index_.first_match(rotated.row(i).transpose(), sorted_.type(i))
        .has_value();
  });
}

bool OverlapChecker::check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot) const {
  Positions const rotated = transformed(sorted_.positions(), rot, trans);
  if (!possible_overlap(rotated)) {
    return false;
  }

  // Greedy bipartite matching, original atoms in order, each taking the
  // lowest-index not-yet-taken image that coincides with it. Collecting the
  // images per original with the image index as the outer loop delivers each
  // original's candidates already ascending.
  auto const n = static_cast<std::size_t>(sorted_.size());
  std::vector<boost::container::small_vector<int, 2>> images(n);
  for (Index ir = 0; ir < sorted_.size(); ++ir) {
    for (int io :
         index_.matches(rotated.row(ir).transpose(), sorted_.type(ir))) {
      images[static_cast<std::size_t>(io)].push_back(static_cast<int>(ir));
    }
  }

  std::vector<bool> taken(n, false);
  return std::ranges::all_of(images, [&](auto const &candidates) {
    auto const free = std::ranges::find_if(candidates, [&](int ir) {
      return !taken[static_cast<std::size_t>(ir)];
    });
    if (free == candidates.end()) {
      return false;
    }
    taken[static_cast<std::size_t>(*free)] = true;
    return true;
  });
}

} // namespace cppcrystal
