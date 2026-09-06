#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/types.hpp>

#include <ranges>
#include <utility>

#pragma GCC visibility push(default)

namespace seitz {

// A crystal cell: a Lattice, fractional positions (row i = atom i), integer
// types, and the per-axis periodicity saying which directions repeat (3D /
// layer / rod / cluster). Immutable: derive with a with_* builder.
class Cell {
public:
  Cell() = default;
  Cell(Lattice lattice, Positions positions, Types types,
       CellPeriodicity periodicity = all_periodic())
      : lattice_{std::move(lattice)}, positions_{std::move(positions)},
        types_{std::move(types)}, periodicity_{periodicity} {}

  [[nodiscard]] Index size() const noexcept {
    return static_cast<Index>(types_.size());
  }

  [[nodiscard]] Lattice const &lattice() const noexcept { return lattice_; }
  [[nodiscard]] Positions const &positions() const noexcept {
    return positions_;
  }
  [[nodiscard]] Types const &types() const noexcept { return types_; }
  [[nodiscard]] CellPeriodicity const &periodicity() const noexcept {
    return periodicity_;
  }

  // Fractional position / type of atom i.
  [[nodiscard]] Vector3d position(Index i) const noexcept {
    return positions_.row(i).transpose();
  }
  [[nodiscard]] int type(Index i) const noexcept {
    return types_[static_cast<std::size_t>(i)];
  }

  // The atoms as {position, type} pairs, so the two containers are never
  // indexed in lockstep at a call site.
  [[nodiscard]] auto atoms() const {
    return std::views::iota(Index{0}, size()) |
           std::views::transform([this](Index i) {
             return std::pair<Vector3d, int>{position(i), type(i)};
           });
  }

  // Same atoms in a new basis, or the same geometry at a new periodicity.
  [[nodiscard]] Cell with_lattice(Lattice lattice) const {
    return Cell{std::move(lattice), positions_, types_, periodicity_};
  }
  [[nodiscard]] Cell with_periodicity(CellPeriodicity periodicity) const {
    return Cell{lattice_, positions_, types_, periodicity};
  }

  // The same crystal in the basis (a' b' c') = (a b c) . basis with its origin
  // at `origin` (fractional, this cell): x' = basis^-1 (x - origin), every atom
  // repeated over the |det basis| lattice points of the new cell, folded along
  // the periodic axes. Errors e_invalid_transformation for a singular matrix
  // or one that mixes an aperiodic axis into a periodic one (an aperiodic
  // axis may only be kept or flipped).
  [[nodiscard]] Result<Cell>
  transformed(Matrix3i const &basis,
              Vector3d const &origin = Vector3d::Zero()) const;

  // The diagonal supercell n[0] x n[1] x n[2].
  [[nodiscard]] Result<Cell> supercell(Vector3i const &n) const {
    return transformed(n.asDiagonal());
  }

  // Every atom shifted by the fractional `shift`, folded along the periodic
  // axes.
  [[nodiscard]] Cell translated(Vector3d const &shift) const;

private:
  Lattice lattice_{};
  Positions positions_{};
  Types types_{};
  CellPeriodicity periodicity_{all_periodic()};
};

} // namespace seitz

#pragma GCC visibility pop
