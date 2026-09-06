#pragma once

#include <seitz/core/cell.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>
#include <seitz/data/element_data.hpp>

#include <boost/container/flat_map.hpp>

#include <utility>

#pragma GCC visibility push(default)

namespace seitz::generate {

// The minimum-distance criterion of a structure, per pair of atom types: an
// explicit override where one was set, otherwise scale * (r(a) + r(b)) over one
// tabulated radius family, an element without that radius taking
// `fallback_radius`. Symmetric in the pair.
class DistanceTolerance {
public:
  // Covalent radii at 0.7 -- the criterion the generators always used.
  DistanceTolerance() = default;

  [[nodiscard]] static DistanceTolerance preset(data::RadiusKind kind,
                                                double scale = 0.7,
                                                double fallback_radius = 1.0) {
    DistanceTolerance out;
    out.kind_ = kind;
    out.scale_ = scale;
    out.fallback_ = fallback_radius;
    return out;
  }

  // Pin the minimum distance of one type pair, whatever the radii say.
  DistanceTolerance &set(int a, int b, double min_distance) {
    overrides_.insert_or_assign(std::minmax(a, b), min_distance);
    return *this;
  }

  // The radius the criterion uses for `type`.
  [[nodiscard]] double radius(int type) const noexcept {
    return data::radius(kind_, type).value_or(fallback_);
  }

  [[nodiscard]] double min_distance(int a, int b) const noexcept {
    if (auto const it = overrides_.find(std::minmax(a, b));
        it != overrides_.end()) {
      return it->second;
    }
    return scale_ * (radius(a) + radius(b));
  }

private:
  data::RadiusKind kind_ = data::RadiusKind::covalent;
  double scale_ = 0.7;
  double fallback_ = 1.0;
  boost::container::flat_map<std::pair<int, int>, double> overrides_;
};

// Smallest Cartesian distance between fractional `a` and `b`: periodic
// components folded to [-0.5, 0.5] with a neighbour search along those axes,
// aperiodic ones kept as the raw difference. `Images::nontrivial` skips the
// zero offset -- for an atom against its own images (a == b).
enum class Images { all, nontrivial };

[[nodiscard]] double minimum_image_distance(
    Vector3d const &a, Vector3d const &b, Matrix3d const &lattice,
    CellPeriodicity const &periodicity, Images images = Images::all) noexcept;

// True iff every atom pair in `cell` (each atom against its own images
// included) clears its type-pair minimum distance, under the cell's
// periodicity.
[[nodiscard]] bool distances_valid(Cell const &cell,
                                   DistanceTolerance const &tol = {});

} // namespace seitz::generate

#pragma GCC visibility pop
