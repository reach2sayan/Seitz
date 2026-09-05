#include "core/overlap.hpp"

#include <algorithm>
#include <cstddef>
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

} // namespace

OverlapChecker::OverlapChecker(Cell const &cell, double symprec)
    : sorted_{sorted_by_distance(cell)}, index_{sorted_, symprec},
      rotated_(sorted_.size(), 3),
      images_(static_cast<std::size_t>(sorted_.size())),
      taken_(static_cast<std::size_t>(sorted_.size()), 0) {}

bool OverlapChecker::possible_overlap(Matrix3d const &rot_transposed,
                                      Vector3d const &trans) const {
  Index const probes = std::min<Index>(sorted_.size(), 3);
  return std::ranges::all_of(std::views::iota(Index{0}, probes), [&](Index i) {
    Vector3d const image =
        (sorted_.positions().row(i) * rot_transposed).transpose() + trans;
    return index_.first_match(image, sorted_.type(i), scratch_).has_value();
  });
}

bool OverlapChecker::check_total_overlap(Vector3d const &trans,
                                         Matrix3i const &rot) const {
  Matrix3d const rot_transposed = rot.cast<double>().transpose();

  // Reject on the probe atoms before mapping the whole cell. Most candidate
  // translations fail here, and building all n images to look at three of them
  // was the single largest cost in the search.
  if (!possible_overlap(rot_transposed, trans)) {
    return false;
  }

  // x -> rot . x + trans applied to every row, into the reused buffer.
  rotated_.noalias() = sorted_.positions() * rot_transposed;
  rotated_.rowwise() += trans.transpose();

  // Greedy bipartite matching, original atoms in order, each taking the
  // lowest-index not-yet-taken image that coincides with it. Collecting the
  // images per original with the image index as the outer loop delivers each
  // original's candidates already ascending.
  for (auto &candidates : images_) {
    candidates.clear();
  }
  for (Index ir = 0; ir < sorted_.size(); ++ir) {
    bool mapped = false;
    for (int io : index_.matches(rotated_.row(ir).transpose(), sorted_.type(ir),
                                 scratch_)) {
      images_[static_cast<std::size_t>(io)].push_back(static_cast<int>(ir));
      mapped = true;
    }
    // An image that coincides with no original can never be claimed, so the n
    // originals would have to be matched from at most n-1 images -- by
    // pigeonhole the matching below is already lost. Bailing here skips the
    // rest of the scan for a translation that was going to fail anyway.
    if (!mapped) {
      return false;
    }
  }

  std::ranges::fill(taken_, std::uint8_t{0});
  return std::ranges::all_of(images_, [&](auto const &candidates) {
    auto const free = std::ranges::find_if(candidates, [&](int ir) {
      return taken_[static_cast<std::size_t>(ir)] == 0;
    });
    if (free == candidates.end()) {
      return false;
    }
    taken_[static_cast<std::size_t>(*free)] = 1;
    return true;
  });
}

} // namespace cppcrystal
