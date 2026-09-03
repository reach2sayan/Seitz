#include <cppcrystal/standardize.hpp>

#include "core/centering.hpp"
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/dataset.hpp>
#include "math/fractional.hpp"
#include "symmetry/primitive.hpp"

#include <utility>

// Cell standardization (3D space-group path). Everything is derived from
// get_dataset: the idealized cases read the dataset's std_* cell directly; the
// no_idealize cases transform the *input* cell into the standardized
// basis/centering (via the dataset transformation matrix), so the input's real
// geometry is preserved.
namespace cppcrystal {

namespace {

using data::Centering;

// Bring `cell` into the primitive setting implied by `trans_mat` (cell ->
// conventional) and the centering (conventional -> primitive), then fold the
// atoms into the primitive lattice.
//   prim_lattice = cell.lattice . trans_mat^-1 . M^-1(centering)
[[nodiscard]] Result<Cell> transform_to_primitive(Cell const &cell,
                                                  Matrix3d const &trans_mat,
                                                  Centering centering,
                                                  double symprec) {
  Lattice const prim_lattice{cell.lattice().matrix() * trans_mat.inverse() *
                             centering_matrix_inv(centering)};
  symmetry::PrimitiveFinder<GroupFamily::space> const finder(
      cell, {symprec, std::nullopt});
  auto trimmed = finder.trim_to(prim_lattice);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

// Expand a primitive cell into the centered conventional cell.
//   conv_lattice  = prim.lattice . M(centering)
//   x_conventional = M^-1(centering) . x_primitive, replicated over the
//   non-trivial centering translations.
[[nodiscard]] Result<Cell> transform_from_primitive(Cell const &primitive,
                                                    Centering centering,
                                                    double symprec) {
  Matrix3d const conv_lattice =
      primitive.lattice().matrix() * centering_matrix(centering).cast<double>();
  Matrix3d const to_conv = centering_matrix_inv(centering);
  auto const shifts = centering_shifts(centering);
  Index const np = primitive.size();
  auto const multi = static_cast<Index>(shifts.size() + 1);

  std::vector<Vector3d> pos;
  Types types;
  pos.reserve(static_cast<std::size_t>(np * multi));
  types.reserve(static_cast<std::size_t>(np * multi));
  for (Index i = 0; i < np; ++i) {
    Vector3d const base = math::wrap_to_unit_cell(
        Vector3d(to_conv * primitive.position(i)));
    pos.push_back(base);
    for (Vector3d const &shift : shifts) {
      pos.push_back(math::wrap_to_unit_cell(Vector3d(base + shift)));
    }
    types.insert(types.end(), static_cast<std::size_t>(multi),
                 primitive.type(i));
  }

  Cell const expanded(Lattice{conv_lattice}, to_positions(pos), std::move(types));
  // Fold/normalize the expanded cell.
  symmetry::PrimitiveFinder<GroupFamily::space> const finder(
      expanded, {symprec, std::nullopt});
  auto trimmed = finder.trim_to(Lattice{conv_lattice});
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

} // namespace

Result<Cell> standardize_cell(Cell const &cell, CellSetting setting,
                              Idealize idealize, Tolerance const &tol) {
  BOOST_LEAF_AUTO(ds, get_dataset(cell, tol));
  Centering const centering = data::spacegroup_type(ds.hall_number).centering;
  double const symprec = tol.symprec;

  if (idealize == Idealize::yes) {
    // The dataset already holds the standardized conventional cell.
    Cell std_cell(Lattice{ds.std_lattice}, ds.std_positions, ds.std_types);
    if (setting == CellSetting::conventional) {
      return std_cell;
    }
    // The bravais cell -> primitive (identity tmat).
    return transform_to_primitive(std_cell, Matrix3d::Identity(), centering,
                                  symprec);
  }

  // Idealize::no: transform the input cell, preserving its real geometry.
  BOOST_LEAF_AUTO(primitive,
                  transform_to_primitive(cell, ds.transformation_matrix,
                                         centering, symprec));
  if (setting == CellSetting::primitive || centering == Centering::primitive) {
    return primitive;
  }
  return transform_from_primitive(primitive, centering, symprec);
}

} // namespace cppcrystal
