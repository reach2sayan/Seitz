#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Everything declared below is the installed ABI: the library is compiled
// with hidden visibility (see CMakeLists.txt), so a public header opens the
// window and closes it again at the end of the file.
#pragma GCC visibility push(default)

namespace cppcrystal::group {

class SpaceGroup;

namespace detail {
struct WyckoffFactory;
}

// One Wyckoff position: a class of points sharing an orbit type under the
// group, parameterised by an affine locus `origin + span(basis columns)` in the
// fractional cell.
//
// The same shape covers both sources. A space- or layer-group position comes
// from the site-symmetry database, where the locus is read off the
// representative coordinate operator and the site-symmetry symbol is tabulated.
// A point- or rod-group position is derived from the arrangement of the
// operations' fixed subspaces, and has no tabulated symbol. Point groups are
// the origin-anchored special case (origin == 0, translation-free operations);
// rod groups carry translations and fold their one periodic axis.
class Wyckoff {
public:
  [[nodiscard]] int multiplicity() const noexcept { return multiplicity_; }

  // Free coordinates of the locus (0..3); the general position has the full
  // periodic dimension.
  [[nodiscard]] int degrees_of_freedom() const noexcept { return dof_; }

  // Wyckoff letter, 'a' = the most special (smallest multiplicity).
  [[nodiscard]] char letter() const noexcept { return letter_; }

  // Tabulated site-symmetry symbol; empty for a derived (point/rod) position,
  // which has no database entry to read one from.
  [[nodiscard]] std::string_view site_symmetry() const noexcept {
    return site_symmetry_;
  }

  // The site-symmetry group: the operations fixing a generic point of the
  // locus. Its order times the multiplicity equals the group's order.
  [[nodiscard]] std::span<SymmetryOperation const> operations() const noexcept {
    return stabiliser_;
  }

  // A point on the locus from `degrees_of_freedom()` free parameters:
  // origin + sum_i params[i] * basis.col(i). Extra parameters are ignored; a
  // 0-DOF position returns its fixed point regardless.
  [[nodiscard]] Vector3d sample(std::span<double const> params) const {
    auto const idx =
        std::views::iota(0, std::min(dof_, static_cast<int>(params.size())));
    return std::ranges::fold_left(
        idx, Vector3d{origin_}, [&](Vector3d acc, int i) -> Vector3d {
          return acc + params[static_cast<std::size_t>(i)] * basis_.col(i);
        });
  }

  // Project a point onto the locus, so a coordinate that is only approximately
  // on the position still yields the exact orbit.
  [[nodiscard]] Vector3d canonical(Vector3d const &xyz) const {
    return projector_ * xyz + projector_shift_;
  }

  // The full orbit of `xyz` (one row per image), projected onto the locus
  // first. Only the periodic axes are folded into the cell: folding an
  // aperiodic one would send a c-flipping image of a layer to 1-z instead of
  // -z. Row count equals multiplicity() for a generic `xyz`.
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

// The shared face of every group family (SpaceGroup, PointGroup, RodGroup):
// a number, a symbol, the operations, and the Wyckoff positions. Plain
// protected inheritance of shared state — there is no runtime hierarchy, no
// virtual dispatch and no polymorphism; the families stay unrelated concrete
// types that happen to answer the same questions.
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
  // The Wyckoff positions, ordered by ascending letter ('a' = the most
  // special); the last is the general position.
  [[nodiscard]] std::span<Wyckoff const> wyckoffs() const noexcept {
    return positions_;
  }

  // Look up a position by letter. Errors if the letter is out of range here.
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

} // namespace cppcrystal::group

#pragma GCC visibility pop
