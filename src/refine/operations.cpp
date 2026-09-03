#include "refine/refinement.hpp"
#include <cppcrystal/core/operation_set.hpp>

#include "core/matrix_order.hpp"
#include "core/position_index.hpp"
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include <cppcrystal/data/spg_database.hpp>

#include <algorithm>
#include <ranges>

// Refined-operation recovery chain (3D path).
namespace cppcrystal::refine {

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

// Bounding-box extent of the parallelepiped spanned by the columns of the
// integer transformation matrix. Each corner is an independent 0/1 choice per
// column, so the per-row max (min) over the corners is the sum of the
// positive (negative) entries of that row.
[[nodiscard]] Vector3i surrounding_frame(Matrix3i const &t_mat) {
  return t_mat.cwiseMax(0).rowwise().sum() - t_mat.cwiseMin(0).rowwise().sum();
}

// Every integer point of the surrounding frame mapped by T^-1 and folded into
// the input cell.
[[nodiscard]] std::vector<Vector3d>
lattice_translations(Vector3i const &frame, Matrix3d const &inv_tmat) {
  std::vector<Vector3d> out;
  out.reserve(static_cast<std::size_t>(frame[0]) *
              static_cast<std::size_t>(frame[1]) *
              static_cast<std::size_t>(frame[2]));

  for (auto const [i, j, k] : std::views::cartesian_product(
           std::views::iota(0, frame[0]), std::views::iota(0, frame[1]),
           std::views::iota(0, frame[2]))) {
    out.emplace_back(
        math::wrap_to_unit_cell(Vector3d(inv_tmat * Vector3d(i, j, k))));
  }

  return out;
}

// Remove overlapping lattice points.
[[nodiscard]] std::vector<Vector3d>
unique_translations(Matrix3d const &lattice,
                    std::vector<Vector3d> const &candidates, double symprec) {
  std::vector<Vector3d> out;
  out.reserve(candidates.size());
  for (const auto &t : candidates) {
    push_unique(out, t, [&](Vector3d const &kept, Vector3d const &cand) {
      return coincident(cand, kept, lattice, symprec, all_periodic());
    });
  }

  return out;
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
  std::vector<Vector3d> const pure_trans = unique_translations(
      lattice, lattice_translations(surrounding_frame(t_mat), inv_tmat),
      symprec);
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
  Operations const conv_sym = operations_from_database(sg.hall);

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

} // namespace cppcrystal::refine
