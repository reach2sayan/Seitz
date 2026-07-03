#include <cppcrystal/refine/operations.hpp>

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>

#include <algorithm>

// Refined-operation recovery chain (3D path).
namespace cppcrystal::refine {

using data::operations_from_database;

namespace {

// t' = t + (R - E) w (the database ops in the shifted origin).
[[nodiscard]] SymmetryOperations
with_origin_shift(SymmetryOperations const &conv_sym, Vector3d const &shift) {
  SymmetryOperations out;
  out.reserve(conv_sym.size());
  std::ranges::transform(
      conv_sym, std::back_inserter(out), [&](const auto &op) {
        const Matrix3d r_minus_e =
            op.rotation.template cast<double>() - Matrix3d::Identity();

        return SymmetryOperation{op.rotation,
                                 Vector3d(op.translation + r_minus_e * shift)};
      });

  return out;
}

// Keep the first operation of each distinct rotation and transform it to the
// primitive setting (R' = T R T^-1, t' = T t).
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

// Bounding-box extent of the parallelepiped spanned by the columns of the
// integer transformation matrix.
[[nodiscard]] Vector3i surrounding_frame(Matrix3i const &t_mat) {
  Eigen::Matrix<int, 3, 8> corners;
  corners.col(0).setZero();
  corners.col(1) = t_mat.col(0);
  corners.col(2) = t_mat.col(1);
  corners.col(3) = t_mat.col(2);
  corners.col(4) = t_mat.col(1) + t_mat.col(2);
  corners.col(5) = t_mat.col(2) + t_mat.col(0);
  corners.col(6) = t_mat.col(0) + t_mat.col(1);
  corners.col(7) = t_mat.col(0) + t_mat.col(1) + t_mat.col(2);

  return corners.rowwise().maxCoeff() - corners.rowwise().minCoeff();
}

// Every integer point of the surrounding frame mapped by T^-1 and folded into
// the input cell.
[[nodiscard]] std::vector<Vector3d>
lattice_translations(Vector3i const &frame, Matrix3d const &inv_tmat) {
  std::vector<Vector3d> out;
  out.reserve(static_cast<std::size_t>(frame[0]) *
              static_cast<std::size_t>(frame[1]) *
              static_cast<std::size_t>(frame[2]));

  for (int i = 0; i < frame[0]; ++i) {
    for (int j = 0; j < frame[1]; ++j) {
      for (int k = 0; k < frame[2]; ++k) {
        out.emplace_back(math::wrap_to_unit_cell(Vector3d(inv_tmat * Vector3d(i, j, k))));
      }
    }
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
    const bool duplicate = std::ranges::any_of(out, [&](const Vector3d &u) {
      return is_overlap(t, u, lattice, symprec);
    });

    if (!duplicate) {
      out.emplace_back(t);
    }
  }

  return out;
}

// Transform the primitive ops to the input cell (R' = T^-1 R T), keeping only
// those that map the input lattice exactly.
[[nodiscard]] SymmetryOperations
symmetry_in_original_cell(Matrix3i const &t_mat, Matrix3d const &inv_tmat,
                          Matrix3d const &lattice,
                          SymmetryOperations const &prim_sym, double symprec) {
  SymmetryOperations out;
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
  return out;
}

// Copy the operations onto every lattice point (fold all axes).
[[nodiscard]] SymmetryOperations
upon_lattice_points(std::vector<Vector3d> const &pure_trans,
                    SymmetryOperations const &t_sym) {
  SymmetryOperations out;
  out.reserve(pure_trans.size() * t_sym.size());
  for (Vector3d const &p : pure_trans) {
    for (auto const &op : t_sym) {
      out.emplace_back(op.rotation, math::wrap_to_unit_cell(Vector3d(op.translation + p)));
    }
  }
  out.shrink_to_fit();
  return out;
}

// Recover the symmetry operations in the original input cell.
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

} // namespace cppcrystal::refine
