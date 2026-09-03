#include "refine/refinement.hpp"

#include "core/centering.hpp"
#include "math/fractional.hpp"
#include "symmetry/primitive.hpp"

#include <utility>
#include <vector>

// The two changes of basis that turn a determined cell into a standardized one
// in a chosen setting. They belong to the refinement because the centering they
// use is the matched group's.
namespace cppcrystal::refine {

template <GroupFamily F>
Result<Cell> Refinement<F>::to_primitive(Cell const &cell,
                                         Matrix3d const &transformation) const {
  //   prim_lattice = cell.lattice . transformation^-1 . M^-1(centering)
  Lattice const prim_lattice{cell.lattice().matrix() *
                             transformation.inverse() *
                             centering_matrix_inv(matched_.type().centering)};
  symmetry::PrimitiveFinder<F> const finder(cell, tol_);
  auto trimmed = finder.trim_to(prim_lattice);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

template <GroupFamily F>
Result<Cell> Refinement<F>::from_primitive(Cell const &primitive) const {
  //   conv_lattice   = prim.lattice . M(centering)
  //   x_conventional = M^-1(centering) . x_primitive, replicated over the
  //   non-trivial centering translations.
  data::Centering const centering = matched_.type().centering;
  Lattice const conv_lattice{primitive.lattice().matrix() *
                             centering_matrix(centering).cast<double>()};
  Matrix3d const to_conv = centering_matrix_inv(centering);
  auto const shifts = centering_shifts(centering);

  std::vector<Vector3d> pos;
  Types types;
  for (auto const &[position, type] : primitive.atoms()) {
    Vector3d const base = math::wrap_to_unit_cell(Vector3d(to_conv * position));
    pos.push_back(base);
    for (Vector3d const &shift : shifts) {
      pos.push_back(math::wrap_to_unit_cell(Vector3d(base + shift)));
    }
    types.insert(types.end(), shifts.size() + 1, type);
  }

  // Fold/normalize the expanded cell.
  Cell const expanded(conv_lattice, to_positions(pos), std::move(types),
                      primitive.periodicity());
  symmetry::PrimitiveFinder<F> const finder(expanded, tol_);
  auto trimmed = finder.trim_to(conv_lattice);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

template Result<Cell>
Refinement<GroupFamily::space>::to_primitive(Cell const &,
                                             Matrix3d const &) const;
template Result<Cell>
Refinement<GroupFamily::layer>::to_primitive(Cell const &,
                                             Matrix3d const &) const;
template Result<Cell>
Refinement<GroupFamily::space>::from_primitive(Cell const &) const;
template Result<Cell>
Refinement<GroupFamily::layer>::from_primitive(Cell const &) const;

} // namespace cppcrystal::refine
