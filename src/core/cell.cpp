#include <seitz/core/cell.hpp>

#include "math/integer_matrix.hpp"

#include <algorithm>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

namespace seitz {

Result<Cell> Cell::transformed(Matrix3i const &basis,
                               Vector3d const &origin) const {
  int const det = math::determinant(math::as_rows(basis));
  // An aperiodic axis has no lattice to change: its row and column must be
  // the (possibly flipped) unit vector.
  auto const untouched = [&](auto axis) {
    auto const i = static_cast<Index>(axis);
    return basis.row(i).cwiseAbs() == Vector3i::Unit(i).transpose() &&
           basis.col(i).cwiseAbs() == Vector3i::Unit(i);
  };
  auto aperiodic = periodicity_ | std::views::enumerate |
                   std::views::filter([](auto const &axis) {
                     return std::get<1>(axis) == AxisKind::aperiodic;
                   }) |
                   std::views::keys;
  if (det == 0 || !std::ranges::all_of(aperiodic, untouched)) {
    return leaf::new_error(e_invalid_transformation{det});
  }

  Matrix3d const inverse = basis.cast<double>().inverse();
  std::vector<Vector3d> const points = math::lattice_points_in_cell(basis);
  auto images = std::views::cartesian_product(atoms(), points);
  std::vector<Vector3d> const rows{
      std::from_range, images | std::views::transform([&](auto const &image) {
        auto const &[atom, point] = image;
        return wrap(Vector3d(inverse * (atom.first - origin) + point),
                    periodicity_);
      })};
  Types types{std::from_range,
              images | std::views::transform([](auto const &image) {
                return std::get<0>(image).second;
              })};
  return Cell{lattice_.transformed(basis.cast<double>()), to_positions(rows),
              std::move(types), periodicity_};
}

Cell Cell::translated(Vector3d const &shift) const {
  std::vector<Vector3d> const rows{
      std::from_range, atoms() | std::views::transform([&](auto const &atom) {
        return wrap(Vector3d(atom.first + shift), periodicity_);
      })};
  return Cell{lattice_, to_positions(rows), types_, periodicity_};
}

} // namespace seitz
