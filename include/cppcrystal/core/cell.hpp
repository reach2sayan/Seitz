#pragma once

#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/types.hpp>

#include <ranges>
#include <utility>

namespace cppcrystal {

// A crystal cell: a Lattice, fractional atomic positions (row i is atom i), and
// integer atom types, together with the per-axis periodicity that says which of
// the three directions actually repeat (3D / layer / rod / cluster).
//
// Immutable: every derived cell is built with a constructor or a with_* builder
// rather than mutated in place, so a Cell can be shared without a copy.
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

  // The atoms as a range of {position, type} pairs, so the two parallel
  // containers are never indexed in lockstep at a call site.
  [[nodiscard]] auto atoms() const {
    return std::views::iota(Index{0}, size()) |
           std::views::transform([this](Index i) {
             return std::pair<Vector3d, int>{position(i), type(i)};
           });
  }

  // Builders for the derived cells the pipeline produces: same atoms in a new
  // basis, or the same geometry reinterpreted at a different periodicity.
  [[nodiscard]] Cell with_lattice(Lattice lattice) const {
    return Cell{std::move(lattice), positions_, types_, periodicity_};
  }
  [[nodiscard]] Cell with_periodicity(CellPeriodicity periodicity) const {
    return Cell{lattice_, positions_, types_, periodicity};
  }

private:
  Lattice lattice_{};
  Positions positions_{};
  Types types_{};
  CellPeriodicity periodicity_{all_periodic()};
};

} // namespace cppcrystal
