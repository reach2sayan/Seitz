#include "refine/refinement.hpp"
#include <seitz/core/operation_set.hpp>

#include "core/matrix_order.hpp"
#include <seitz/core/fractional.hpp>
#include "math/integer_matrix.hpp"
#include <seitz/data/spg_database.hpp>

#include <ranges>

// Refined-operation recovery chain (3D path).
namespace seitz::refine {

using data::operations_from_database;

namespace {

// t' = t + (R - E) w (the database ops in the shifted origin).
[[nodiscard]] Operations with_origin_shift(Operations const &conv_sym,
                                           Vector3d const &shift) {
  return Operations{
      std::from_range,
      conv_sym | std::views::transform([&](auto const &op) {
        Matrix3d const r_minus_e =
            op.rotation.template cast<double>() - Matrix3d::Identity();
        return SymmetryOperation{
            .rotation = op.rotation,
            .translation = Vector3d(op.translation + r_minus_e * shift)};
      })};
}

// Keep the first operation of each distinct rotation and transform it to the
// primitive setting (R' = T R T^-1, t' = T t).
[[nodiscard]] Operations primitive_db_symmetry(Matrix3d const &t_mat,
                                               Operations const &conv_sym) {
  return Operations{unique_by_rotation(conv_sym, &SymmetryOperation::rotation)}
      .conjugated_by(t_mat, t_mat.inverse());
}

// Transform the primitive ops to the input cell (R' = T^-1 R T), keeping only
// those that map the input lattice exactly.
[[nodiscard]] Operations symmetry_in_original_cell(Matrix3i const &t_mat,
                                                   Matrix3d const &inv_tmat,
                                                   Matrix3d const &lattice,
                                                   Operations const &prim_sym,
                                                   double symprec) {
  std::vector<SymmetryOperation> out;
  out.reserve(prim_sym.size());
  for (auto const &op : prim_sym) {
    Matrix3d const rot_d =
        inv_tmat * op.rotation.cast<double>() * t_mat.cast<double>();
    Matrix3i const rot_i = math::round_to_int(rot_d);
    Matrix3d const lat_i = lattice * rot_i.cast<double>();
    Matrix3d const lat_d = lattice * rot_d;
    if ((lat_i - lat_d).cwiseAbs().maxCoeff() <= symprec) {
      out.emplace_back(rot_i, Vector3d(inv_tmat * op.translation));
    }
  }
  out.shrink_to_fit();
  return Operations{std::move(out)};
}

// Copy the operations onto every lattice point (fold all axes).
[[nodiscard]] Operations
upon_lattice_points(std::vector<Vector3d> const &pure_trans,
                    Operations const &t_sym) {
  std::vector<SymmetryOperation> out;
  out.reserve(pure_trans.size() * t_sym.size());
  for (Vector3d const &p : pure_trans) {
    for (auto const &op : t_sym) {
      out.emplace_back(op.rotation,
                       math::wrap_to_unit_cell(Vector3d(op.translation + p)));
    }
  }
  out.shrink_to_fit();
  return Operations{std::move(out)};
}

// Recover the symmetry operations in the original input cell.
[[nodiscard]] std::optional<Operations>
recover_in_original_cell(Operations const &prim_sym, Matrix3i const &t_mat,
                         Matrix3d const &lattice, int multiplicity,
                         double symprec) {
  Matrix3d const inv_tmat = t_mat.cast<double>().inverse();
  std::vector const pure_trans = math::lattice_points_in_cell(t_mat);
  Operations const t_sym =
      symmetry_in_original_cell(t_mat, inv_tmat, lattice, prim_sym, symprec);

  if (static_cast<int>(pure_trans.size()) != multiplicity) {
    return std::nullopt;
  }
  return upon_lattice_points(pure_trans, t_sym);
}

} // namespace

template <GroupFamily F>
std::optional<Operations> Refinement<F>::operations() const {
  SpacegroupMatch const &sg = matched_;
  Cell const &primitive = primitive_;
  Cell const &cell = cell_;
  double const symprec = tol_.symprec;
  Operations const &conv_sym = operations_from_database(sg.hall);

  Matrix3d const inv_prim = primitive.lattice().matrix().inverse();

  // Conventional ops in the shifted origin, transformed to the primitive cell.
  Matrix3d const t_mat_pb = inv_prim * sg.bravais_lattice;
  Operations const prim_sym = primitive_db_symmetry(
      t_mat_pb, with_origin_shift(conv_sym, sg.origin_shift));

  // Recover them in the input cell.
  Matrix3i const t_mat_pc =
      math::round_to_int(Matrix3d(inv_prim * cell.lattice().matrix()));
  int const multiplicity = static_cast<int>(cell.size() / primitive.size());
  return recover_in_original_cell(prim_sym, t_mat_pc, cell.lattice().matrix(),
                                  multiplicity, symprec);
}

template std::optional<Operations>
Refinement<GroupFamily::space>::operations() const;
template std::optional<Operations>
Refinement<GroupFamily::layer>::operations() const;

} // namespace seitz::refine
