#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/group/locus_wyckoff.hpp>

#include <algorithm>
#include <concepts>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

// The shared Wyckoff-derivation driver behind the derived group families:
// build the arrangement of stabiliser loci (the whole space, every
// operation's fixed loci, the closure under intersection), emit one derived
// position per orbit of loci under the group, then sort and assign letters.
// A geometry policy supplies everything dimension-specific.
namespace cppcrystal::group::detail {

// One derived Wyckoff position before letter assignment.
struct DerivedLocus {
  int multiplicity = 0;
  int dof = 0;
  Vector3d origin{Vector3d::Zero()};
  Matrix3d locus_basis{Matrix3d::Zero()};
  SymmetryOperations orbit_ops;
  SymmetryOperations site_symmetry;
};

struct WyckoffFactory {
  // Sort ('a' = most special: smallest multiplicity, then fewest degrees of
  // freedom; the general position lands last) and assign letters.
  [[nodiscard]] static std::vector<LocusWyckoff>
  to_wyckoffs(std::vector<DerivedLocus> derived) {
    std::ranges::sort(derived, {}, [](DerivedLocus const &d) {
      return std::pair{d.multiplicity, d.dof};
    });
    std::vector<LocusWyckoff> out;
    out.reserve(derived.size());
    char letter = 'a';
    for (auto &d : derived) {
      out.push_back(LocusWyckoff{d.multiplicity, d.dof, letter++,
                                 std::move(d.origin), std::move(d.locus_basis),
                                 std::move(d.orbit_ops),
                                 std::move(d.site_symmetry)});
    }
    return out;
  }
};

// Orbit/stabiliser split of a group at one point: walk the operations, keep
// the stabilizer of p0 and one coset representative per distinct image.
// `fold` normalises an image into the fundamental domain (wrap_to_unit_cell,
// a single-axis fold, or identity); `same_point` is the geometry's equality.
struct OrbitPartition {
  SymmetryOperations orbit_ops;     // coset representatives, one per point
  SymmetryOperations site_symmetry; // stabilizer of p0
  std::vector<Vector3d> points;     // the distinct (folded) orbit points
};

template <class Fold, class SamePoint>
[[nodiscard]] OrbitPartition
partition_orbit(std::span<SymmetryOperation const> ops, Vector3d const &p0,
                Fold &&fold, SamePoint &&same_point) {
  OrbitPartition out;
  for (auto const &op : ops) {
    Vector3d const image = fold(op.apply(p0));
    if (same_point(image, p0)) {
      out.site_symmetry.push_back(op);
    }
    if (std::ranges::none_of(out.points, [&](Vector3d const &q) {
          return same_point(q, image);
        })) {
      out.points.push_back(image);
      out.orbit_ops.push_back(op);
    }
  }
  return out;
}

// A dimension policy for the locus arrangement. `Locus` is the geometry's
// stabiliser-locus representation; `derive` computes the orbit data of one
// representative locus (the geometry holds the group's operations).
template <class G>
concept LocusGeometry =
    requires(G const g, typename G::Locus const l, SymmetryOperation const &op) {
      { g.whole_space() } -> std::same_as<typename G::Locus>;
      { g.fixed_loci(op) } -> std::ranges::input_range;
      { g.intersections(l, l) } -> std::ranges::input_range;
      { g.same_locus(l, l) } -> std::convertible_to<bool>;
      { g.image(op, l) } -> std::same_as<typename G::Locus>;
      { g.derive(l) } -> std::same_as<DerivedLocus>;
    };

template <LocusGeometry G>
[[nodiscard]] std::vector<LocusWyckoff>
derive_wyckoff_positions(std::span<SymmetryOperation const> ops, G const &geom) {
  using Locus = typename G::Locus;

  // 1. The arrangement: whole space (general position), every operation's
  //    fixed loci, and the closure under pairwise intersection (re-scanning
  //    newly added loci).
  std::vector<Locus> loci;
  auto add_unique = [&](Locus const &candidate) {
    if (std::ranges::none_of(loci, [&](Locus const &l) {
          return geom.same_locus(l, candidate);
        })) {
      loci.push_back(candidate);
    }
  };
  add_unique(geom.whole_space());
  for (auto const &op : ops) {
    for (auto const &l : geom.fixed_loci(op)) {
      add_unique(l);
    }
  }
  for (std::size_t i = 0; i < loci.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      for (auto const &l : geom.intersections(loci[i], loci[j])) {
        add_unique(l);
      }
    }
  }

  // 2. One derived position per orbit of loci (conjugate stabilisers): derive
  //    the representative, then mark every image of its locus as assigned.
  std::vector<DerivedLocus> derived;
  std::vector<bool> assigned(loci.size(), false);
  for (std::size_t i = 0; i < loci.size(); ++i) {
    if (assigned[i]) {
      continue;
    }
    derived.push_back(geom.derive(loci[i]));
    for (auto const &op : ops) {
      Locus const mapped = geom.image(op, loci[i]);
      for (std::size_t k = 0; k < loci.size(); ++k) {
        if (!assigned[k] && geom.same_locus(loci[k], mapped)) {
          assigned[k] = true;
        }
      }
    }
  }

  return WyckoffFactory::to_wyckoffs(std::move(derived));
}

} // namespace cppcrystal::group::detail
