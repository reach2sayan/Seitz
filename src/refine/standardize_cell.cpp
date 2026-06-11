#include <spglib/standardize.hpp>

#include <spglib/core/centering.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/dataset.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/symmetry/primitive.hpp>

#include <utility>

// Port of spglib.c's standardize_cell / standardize_primitive /
// get_standardized_cell (3D space-group path). Everything is derived from
// get_dataset: the idealized cases read the dataset's std_* cell directly; the
// no_idealize cases transform the *input* cell into the standardized
// basis/centering (via the dataset transformation matrix), so the input's real
// geometry is preserved.
namespace spglib {

namespace {

using data::Centering;

// spa_transform_to_primitive: bring `cell` into the primitive setting implied by
// `trans_mat` (cell -> conventional) and the centering (conventional ->
// primitive), then fold the atoms into the primitive lattice.
//   prim_lattice = cell.lattice . trans_mat^-1 . M^-1(centering)
[[nodiscard]] Result<Cell> transform_to_primitive(Cell const &cell,
                                                  Matrix3d const &trans_mat,
                                                  Centering centering,
                                                  double symprec) {
  Matrix3d const prim_lattice =
      cell.lattice() * trans_mat.inverse() * centering_matrix_inv(centering);
  auto trimmed = symmetry::trim_to_lattice(prim_lattice, cell, symprec);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

// spa_transform_from_primitive: expand a primitive cell into the centered
// conventional cell.
//   conv_lattice  = prim.lattice . M(centering)
//   x_conventional = M^-1(centering) . x_primitive, replicated over the
//   non-trivial centering translations.
[[nodiscard]] Result<Cell> transform_from_primitive(Cell const &primitive,
                                                    Centering centering,
                                                    double symprec) {
  Matrix3d const conv_lattice =
      primitive.lattice() * centering_matrix(centering).cast<double>();
  Matrix3d const to_conv = centering_matrix_inv(centering);
  std::vector<Vector3d> const shifts = centering_shifts(centering);
  Index const np = primitive.size();
  auto const multi = static_cast<Index>(shifts.size() + 1);

  Positions pos(np * multi, 3);
  Types types(static_cast<std::size_t>(np * multi));
  Index row = 0;
  for (Index i = 0; i < np; ++i) {
    Vector3d const base = math::wrap_to_unit_cell(
        Vector3d(to_conv * primitive.position(i)));
    pos.row(row) = base.transpose();
    types[static_cast<std::size_t>(row)] = primitive.type(i);
    ++row;
    for (Vector3d const &shift : shifts) {
      pos.row(row) = math::wrap_to_unit_cell(Vector3d(base + shift)).transpose();
      types[static_cast<std::size_t>(row)] = primitive.type(i);
      ++row;
    }
  }

  Cell const expanded(conv_lattice, std::move(pos), std::move(types));
  // Fold/normalize (the reference runs cel_trim_cell on the expanded cell).
  auto trimmed = symmetry::trim_to_lattice(conv_lattice, expanded, symprec);
  if (!trimmed) {
    return leaf::new_error(e_cell_standardization_failed{});
  }
  return std::move(trimmed->first);
}

} // namespace

Result<Cell> standardize_cell(Cell const &cell, StandardizeOptions options,
                              double symprec, AngleTolerance angle_tolerance) {
  BOOST_LEAF_AUTO(ds, get_dataset(cell, symprec, angle_tolerance));
  Centering const centering =
      data::spacegroup_type(ds.hall_number).centering;

  if (!options.no_idealize) {
    // Idealized: the dataset already holds the standardized conventional cell.
    Cell std_cell(ds.std_lattice, ds.std_positions, ds.std_types);
    if (!options.to_primitive) {
      return std_cell;
    }
    // standardize_primitive: the bravais cell -> primitive (identity tmat).
    return transform_to_primitive(std_cell, Matrix3d::Identity(), centering,
                                  symprec);
  }

  // no_idealize: transform the input cell, preserving its real geometry.
  BOOST_LEAF_AUTO(primitive,
                  transform_to_primitive(cell, ds.transformation_matrix,
                                         centering, symprec));
  if (options.to_primitive || centering == Centering::primitive) {
    return primitive;
  }
  return transform_from_primitive(primitive, centering, symprec);
}

} // namespace spglib
