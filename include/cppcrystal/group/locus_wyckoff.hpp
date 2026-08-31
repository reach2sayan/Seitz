#pragma once

#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace cppcrystal::group {

namespace detail {
struct WyckoffFactory;
}

// A Wyckoff position derived from a locus arrangement: a class of points
// sharing an orbit type under the group, parameterised by an affine locus
// `origin + span(basis columns)` in the fractional cell. Point groups are the
// origin-anchored special case (origin == 0, operations translation-free);
// rod groups carry translations and fold their one periodic axis. `sample()`
// turns free parameters into a point on the locus and `orbit_operations()`
// expands that point into the full orbit.
class LocusWyckoff {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Free coordinates of the locus (0..3); the general position has the full
  // periodic dimension.
  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  // Wyckoff letter, 'a' = the most special (smallest multiplicity).
  [[nodiscard]] char letter() const noexcept { return letter_; }

  // The site-symmetry group: the operations fixing a generic point of the
  // locus. Its order times the multiplicity equals the group's order.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return site_symmetry_;
  }

  // Coset representatives generating the orbit (one per orbit point).
  // Applying each to a point from sample() yields the full orbit.
  [[nodiscard]] std::span<SymmetryOperation const>
  orbit_operations() const noexcept {
    return orbit_ops_;
  }

  // A point on the locus (fractional coordinate) from `degrees_of_freedom()`
  // free parameters: origin + sum_i params[i] * basis.col(i). Extra parameters
  // are ignored; a 0-DOF position returns the fixed point regardless.
  [[nodiscard]] Vector3d sample(std::span<double const> params) const {
    auto const idx =
        std::views::iota(0, std::min(dof_, static_cast<int>(params.size())));
    return std::ranges::fold_left(
        idx, Vector3d{locus_origin_}, [&](Vector3d acc, int i) -> Vector3d {
          return acc + params[static_cast<std::size_t>(i)] *
                           locus_basis_.col(i);
        });
  }

  // The i-th free direction of the locus (fractional), 0 <= i < dof. Rod
  // bases are axis-separated: a direction is either the periodic axis or lies
  // in the aperiodic cross-section — generators use this to decide how each
  // free coordinate is sampled.
  [[nodiscard]] Vector3d free_direction(int i) const {
    return locus_basis_.col(i);
  }

private:
  friend struct detail::WyckoffFactory;
  LocusWyckoff(int multiplicity, int dof, char letter, Vector3d locus_origin,
               Matrix3d locus_basis, SymmetryOperations orbit_ops,
               SymmetryOperations site_symmetry)
      : multiplicity_(multiplicity), dof_(dof), letter_(letter),
        locus_origin_(std::move(locus_origin)),
        locus_basis_(std::move(locus_basis)), orbit_ops_(std::move(orbit_ops)),
        site_symmetry_(std::move(site_symmetry)) {}

  int multiplicity_ = 0;
  int dof_ = 0;
  char letter_ = 'a';
  Vector3d locus_origin_{Vector3d::Zero()};
  // Columns 0..dof_-1 are a basis of the locus directions; the rest unused.
  Matrix3d locus_basis_{Matrix3d::Zero()};
  SymmetryOperations orbit_ops_;
  SymmetryOperations site_symmetry_;
};

// Shared state and accessors of the derived group families (PointGroup,
// RodGroup). Plain protected inheritance — there is no runtime hierarchy and
// no virtual dispatch; the families remain unrelated concrete types.
class GroupBase {
public:
  [[nodiscard]] int number() const noexcept { return number_; }
  [[nodiscard]] std::string_view symbol() const noexcept { return symbol_; }

  // The order of the group.
  [[nodiscard]] int order() const noexcept {
    return static_cast<int>(operations_.size());
  }

  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return operations_;
  }

  // The Wyckoff positions, ordered by ascending multiplicity ('a' = the most
  // special); the last is the general position.
  [[nodiscard]] std::span<LocusWyckoff const> wyckoffs() const noexcept {
    return positions_;
  }

protected:
  GroupBase() = default;

  int number_ = 0;
  std::string_view symbol_;
  SymmetryOperations operations_;
  std::vector<LocusWyckoff> positions_;
};

} // namespace cppcrystal::group
