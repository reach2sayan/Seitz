#include <spglib/refine/operations.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/data/spg_database.hpp>
#include <spglib/math/fractional.hpp>
#include <spglib/math/integer_matrix.hpp>

#include <algorithm>
#include <array>

// Port of the refined-operation recovery chain in refinement.c (3D path).
namespace spglib::refine {

using data::operations_from_database;

namespace {

// set_translation_with_origin_shift: t' = t + (R - E) w (the database ops in
// the shifted origin).
[[nodiscard]] SymmetryOperations
with_origin_shift(SymmetryOperations const &conv_sym, Vector3d const &shift) {
  SymmetryOperations out;
  out.reserve(conv_sym.size());
  for (auto const &op : conv_sym) {
    Matrix3d const r_minus_e =
        op.rotation.cast<double>() - Matrix3d::Identity();
    out.push_back({op.rotation, Vector3d(op.translation + r_minus_e * shift)});
  }
  return out;
}

// get_primitive_db_symmetry: keep the first operation of each distinct rotation
// and transform it to the primitive setting (R' = T R T^-1, t' = T t).
[[nodiscard]] SymmetryOperations
primitive_db_symmetry(Matrix3d const &t_mat,
                      SymmetryOperations const &conv_sym) {
  Matrix3d const inv_mat = t_mat.inverse();
  SymmetryOperations out;
  std::vector<Matrix3i> seen_rotations; // first occurrence
  for (auto const &op : conv_sym) {
    if (std::ranges::any_of(seen_rotations, [&](Matrix3i const &r) {
          return r == op.rotation;
        })) {
      continue;
    }
    seen_rotations.push_back(op.rotation);
    Matrix3d const rot_d = t_mat * op.rotation.cast<double>() * inv_mat;
    out.push_back(
        {math::round_to_int(rot_d), Vector3d(t_mat * op.translation)});
  }
  return out;
}

// get_surrounding_frame: bounding-box extent of the parallelepiped spanned by
// the columns of the integer transformation matrix.
[[nodiscard]] Vector3i surrounding_frame(Matrix3i const &t_mat) {
  std::array<Vector3i, 8> const corners = {
      Vector3i::Zero(),
      t_mat.col(0),
      t_mat.col(1),
      t_mat.col(2),
      Vector3i(t_mat.col(1) + t_mat.col(2)),
      Vector3i(t_mat.col(2) + t_mat.col(0)),
      Vector3i(t_mat.col(0) + t_mat.col(1)),
      Vector3i(t_mat.col(0) + t_mat.col(1) + t_mat.col(2))};

  Vector3i frame;
  for (int i = 0; i < 3; ++i) {
    auto const coord = [i](Vector3i const &c) { return c[i]; };
    auto const [min, max] = std::ranges::minmax(corners, {}, coord);
    frame[i] = coord(max) - coord(min);
  }
  return frame;
}

// get_lattice_translations: every integer point of the surrounding frame mapped
// by T^-1 and folded into the input cell.
[[nodiscard]] std::vector<Vector3d>
lattice_translations(Vector3i const &frame, Matrix3d const &inv_tmat) {
  std::vector<Vector3d> out;
  for (int i = 0; i < frame[0]; ++i)
    for (int j = 0; j < frame[1]; ++j)
      for (int k = 0; k < frame[2]; ++k) {
        out.push_back(math::mod1(Vector3d(inv_tmat * Vector3d(i, j, k))));
      }
  return out;
}

// remove_overlapping_lattice_points.
[[nodiscard]] std::vector<Vector3d>
unique_translations(Matrix3d const &lattice,
                    std::vector<Vector3d> const &candidates, double symprec) {
  std::vector<Vector3d> out;
  for (auto const &t : candidates) {
    if (std::ranges::none_of(out, [&](Vector3d const &u) {
          return is_overlap(t, u, lattice, symprec);
        })) {
      out.push_back(t);
    }
  }
  return out;
}

// get_symmetry_in_original_cell: transform the primitive ops to the input cell
// (R' = T^-1 R T), keeping only those that map the input lattice exactly.
[[nodiscard]] SymmetryOperations
symmetry_in_original_cell(Matrix3i const &t_mat, Matrix3d const &inv_tmat,
                          Matrix3d const &lattice,
                          SymmetryOperations const &prim_sym, double symprec) {
  SymmetryOperations out;
  for (auto const &op : prim_sym) {
    Matrix3d const rot_d =
        inv_tmat * op.rotation.cast<double>() * t_mat.cast<double>();
    Matrix3i const rot_i = math::round_to_int(rot_d);
    Matrix3d const lat_i = lattice * rot_i.cast<double>();
    Matrix3d const lat_d = lattice * rot_d;
    if ((lat_i - lat_d).cwiseAbs().maxCoeff() <= symprec) {
      out.push_back({rot_i, Vector3d(inv_tmat * op.translation)});
    }
  }
  return out;
}

// copy_symmetry_upon_lattice_points (aperiodic_axis = -1: fold all axes).
[[nodiscard]] SymmetryOperations
upon_lattice_points(std::vector<Vector3d> const &pure_trans,
                    SymmetryOperations const &t_sym) {
  SymmetryOperations out;
  out.reserve(pure_trans.size() * t_sym.size());
  for (Vector3d const &p : pure_trans) {
    for (auto const &op : t_sym) {
      out.push_back({op.rotation, math::mod1(Vector3d(op.translation + p))});
    }
  }
  return out;
}

// recover_symmetry_in_original_cell.
[[nodiscard]] std::optional<SymmetryOperations>
recover_in_original_cell(SymmetryOperations const &prim_sym,
                         Matrix3i const &t_mat, Matrix3d const &lattice,
                         int multiplicity, double symprec) {
  Matrix3d const inv_tmat = t_mat.cast<double>().inverse();
  std::vector<Vector3d> const pure_trans = unique_translations(
      lattice, lattice_translations(surrounding_frame(t_mat), inv_tmat),
      symprec);
  SymmetryOperations const t_sym =
      symmetry_in_original_cell(t_mat, inv_tmat, lattice, prim_sym, symprec);

  if (static_cast<int>(pure_trans.size()) != multiplicity) {
    return std::nullopt;
  }
  return upon_lattice_points(pure_trans, t_sym);
}

} // namespace

std::optional<SymmetryOperations>
refined_operations(spacegroup::Spacegroup const &sg, Cell const &primitive,
                   Cell const &cell, double symprec) {
  SymmetryOperations const conv_sym =
      operations_from_database(sg.type.hall_number);

  Matrix3d const inv_prim = primitive.lattice().inverse();

  // Conventional ops in the shifted origin, transformed to the primitive cell.
  Matrix3d const t_mat_pb = inv_prim * sg.bravais_lattice;
  SymmetryOperations const prim_sym = primitive_db_symmetry(
      t_mat_pb, with_origin_shift(conv_sym, sg.origin_shift));

  // Recover them in the input cell.
  Matrix3i const t_mat_pc =
      math::round_to_int(Matrix3d(inv_prim * cell.lattice()));
  int const multiplicity = static_cast<int>(cell.size() / primitive.size());
  return recover_in_original_cell(prim_sym, t_mat_pc, cell.lattice(),
                                  multiplicity, symprec);
}

} // namespace spglib::refine
