#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::group {

class SpaceGroup;

namespace detail {
struct WyckoffFactory;
}

// One Wyckoff position: points sharing an orbit type under the group, as the
// affine locus origin + span(basis columns) in the fractional cell.
//
// Space/layer positions come from the site-symmetry database (locus read off
// the representative coordinate operator, symbol tabulated); point/rod ones are
// derived from the arrangement of the operations' fixed subspaces and carry no
// symbol. Point groups are the origin-anchored case (origin = 0,
// translation-free ops); rod groups translate and fold their periodic axis.
class Wyckoff {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Free coordinates of the locus (0..3); the general position has the full
  // periodic dimension.
  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  // Wyckoff letter, 'a' = the most special (smallest multiplicity).
  [[nodiscard]] char letter() const noexcept { return letter_; }

  // Tabulated site-symmetry symbol; empty for a derived (point/rod) position.
  [[nodiscard]] std::string_view site_symmetry() const noexcept {
    return site_symmetry_;
  }

  // Site-symmetry group: the ops fixing a generic point of the locus.
  // |stabiliser| * multiplicity = |G|.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return stabiliser_;
  }

  // origin + sum_i params[i] * basis.col(i) over the free parameters. Extra
  // parameters are ignored; a 0-DOF position returns its fixed point.
  [[nodiscard]] Vector3d sample(std::span<double const> params) const {
    auto const idx =
        std::views::iota(0, std::min(dof_, static_cast<int>(params.size())));
    return std::ranges::fold_left(
        idx, Vector3d{origin_}, [&](Vector3d acc, int i) -> Vector3d {
          return acc + params[static_cast<std::size_t>(i)] * basis_.col(i);
        });
  }

  // Projection onto the locus: an approximate coordinate still gives the
  // exact orbit.
  [[nodiscard]] Vector3d canonical(Vector3d const &xyz) const {
    return projector_ * xyz + projector_shift_;
  }

  // Full orbit of `xyz` (a row per image), projected onto the locus first.
  // Only periodic axes fold into the cell -- folding an aperiodic one sends a
  // c-flipping image of a layer to 1-z instead of -z. Rows = multiplicity()
  // for generic `xyz`.
  [[nodiscard]] Positions
  orbit(Vector3d const &xyz,
        CellPeriodicity const &periodicity = all_periodic()) const;

private:
  friend struct detail::WyckoffFactory;
  friend class SpaceGroup;
  Wyckoff(int multiplicity, int dof, char letter, std::string_view symbol,
          Vector3d origin, Matrix3d basis, Matrix3d projector,
          Vector3d projector_shift, Operations orbit_ops, Operations stabiliser)
      : multiplicity_(multiplicity), dof_(dof), letter_(letter),
        site_symmetry_(symbol), origin_(std::move(origin)),
        basis_(std::move(basis)), projector_(std::move(projector)),
        projector_shift_(std::move(projector_shift)),
        orbit_ops_(std::move(orbit_ops)), stabiliser_(std::move(stabiliser)) {}

  int multiplicity_ = 0;
  int dof_ = 0;
  char letter_ = 'a';
  std::string_view site_symmetry_;
  Vector3d origin_{Vector3d::Zero()};
  // Columns 0..dof_-1 are a basis of the locus directions; the rest unused.
  Matrix3d basis_{Matrix3d::Zero()};
  // canonical(x) = projector_ x + projector_shift_, idempotent on the locus.
  Matrix3d projector_{Matrix3d::Zero()};
  Vector3d projector_shift_{Vector3d::Zero()};
  Operations orbit_ops_;
  Operations stabiliser_;
};

// The shared face of every family (SpaceGroup, PointGroup, RodGroup): number,
// symbol, operations, Wyckoff positions. Protected inheritance of shared state
// only -- no virtuals, no runtime hierarchy; the families stay unrelated
// concrete types answering the same questions.
class GroupBase {
public:
  [[nodiscard]] int number() const noexcept { return number_; }
  [[nodiscard]] std::string_view symbol() const noexcept { return symbol_; }
  [[nodiscard]] int order() const noexcept {
    return static_cast<int>(operations_.size());
  }
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return operations_;
  }
  // Ordered by ascending letter ('a' = most special); the last is general.
  [[nodiscard]] std::span<Wyckoff const> wyckoffs() const noexcept {
    return positions_;
  }

  // Errors if the letter names no position of this group.
  [[nodiscard]] Result<Wyckoff const *> wyckoff(char letter) const {
    if (auto const it = std::ranges::find(positions_, letter, &Wyckoff::letter);
        it != positions_.end()) {
      return &*it;
    }
    return leaf::new_error(
        e_message{std::string("wyckoff: no Wyckoff position '") + letter +
                  "' in this group"});
  }

protected:
  GroupBase() = default;
  int number_ = 0;
  std::string_view symbol_;
  Operations operations_;
  std::vector<Wyckoff> positions_;
};

} // namespace seitz::group

#pragma GCC visibility pop
