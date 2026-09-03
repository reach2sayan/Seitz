#include <cppcrystal/kpoint/reciprocal_mesh.hpp>

#include "core/matrix_order.hpp" // unique_by_rotation
#include <cppcrystal/core/tolerance.hpp>    // approx_equal
#include <cppcrystal/analysis/symmetry_analyzer.hpp>           // get_dataset
#include <cppcrystal/kpoint/grid.hpp>
#include "math/fractional.hpp" // math::nint

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <vector>

namespace cppcrystal::kpoint {

std::vector<Matrix3i>
point_group_reciprocal(std::vector<Matrix3i> const &rotations,
                       TimeReversal time_reversal) {
  // All transposes first, then (with time reversal) all inversion partners
  // (-transpose); de-duplication keeps first occurrence, so order matters.
  std::vector<Matrix3i> candidates;
  candidates.reserve(rotations.size() *
                     (time_reversal == TimeReversal::on ? 2 : 1));

  std::ranges::transform(rotations, std::back_inserter(candidates),
                         [](const Matrix3i &r) { return r.transpose(); });
  if (time_reversal == TimeReversal::on) {
    std::ranges::transform(rotations, std::back_inserter(candidates),
                           [](const Matrix3i &r) { return -r.transpose(); });
  }

  return unique_by_rotation(candidates);
}

namespace {

// The subset of `rot_reciprocal` that maps the q-point set onto itself (modulo
// a reciprocal lattice vector) within `symprec` — the stabilizer of the
// q-points. Only stabilized_reciprocal_mesh needs it.
[[nodiscard]] std::vector<Matrix3i>
point_group_reciprocal_with_q(std::vector<Matrix3i> const &rot_reciprocal,
                              double symprec,
                              std::vector<Vector3d> const &qpoints) {

  auto preserves_qpoints = [&](const Matrix3i &rot) {
    return std::ranges::all_of(qpoints, [&](const Vector3d &q) {
      const Vector3d q_rot = rot.cast<double>() * q;
      return std::ranges::any_of(qpoints, [&](const Vector3d &other) {
        Vector3d const diff = q_rot - other;
        return approx_equal(diff, math::round_to_int(diff).cast<double>(),
                            symprec);
      });
    });
  };

  std::vector<Matrix3i> out;
  std::ranges::copy(rot_reciprocal | std::views::filter(preserves_qpoints),
                    std::back_inserter(out));
  return out;
}

// True -> the fast "normal" reduction; false -> the "distortion" path (3/6-fold
// rotations or non-conventional cells). NOTE the deliberate quirk that the a=b
// and b=c flags test the same column — kept for bit-identical behaviour against
// the reference.
[[nodiscard]] bool
mesh_has_conventional_symmetry(Vector3i const &mesh, Vector3i const &is_shift,
                               std::vector<Matrix3i> const &rot_reciprocal) {
  if (std::ranges::any_of(rot_reciprocal, [](Matrix3i const &rot) {
        return rot.cwiseAbs().sum() > 3;
      })) {
    return false;
  }

  Vector3i const e_b(0, 1, 0);
  Vector3i const e_c(0, 0, 1);
  bool const eq_ab = std::ranges::any_of( // a == b axis present
      rot_reciprocal, [&](Matrix3i const &rot) { return rot.col(0) == e_b; });
  bool const eq_bc = eq_ab; // b == c (intentionally the same column as a == b)
  bool const eq_ca = std::ranges::any_of( // c == a
      rot_reciprocal, [&](Matrix3i const &rot) { return rot.col(0) == e_c; });

  return (!eq_ab || (mesh[0] == mesh[1] && is_shift[0] == is_shift[1])) &&
         (!eq_bc || (mesh[1] == mesh[2] && is_shift[1] == is_shift[2])) &&
         (!eq_ca || (mesh[2] == mesh[0] && is_shift[2] == is_shift[0]));
}

// Smallest grid-point index in the orbit of `address_double` under the
// reciprocal group (the irreducible representative). The group is closed, so
// iterating every rotation reaches the whole orbit; the minimum is canonical.
[[nodiscard]] std::size_t
normal_representative(Vector3i const &address_double, Vector3i const &mesh,
                      std::size_t self,
                      std::vector<Matrix3i> const &rot_reciprocal) {
  auto const rotated =
      rot_reciprocal | std::views::transform([&](Matrix3i const &rot) {
        return grid_point_double_mesh(Vector3i(rot * address_double), mesh);
      });
  std::size_t best = self;
  for (std::size_t const rot_point : rotated) {
    best = std::min(best, rot_point);
  }
  return best;
}

// As above, but for the distortion path: the rotation is applied in 64-bit
// (scaled by `divisor`) and a rotated point is only valid when it divides
// evenly and its parity matches `is_shift`.
[[nodiscard]] std::size_t distortion_representative(
    Vector3i const &address_double, Vector3i const &mesh,
    Vector3i const &is_shift, std::size_t self,
    std::array<std::int64_t, 3> const &divisor,
    std::vector<Matrix3i> const &rot_reciprocal) {
  std::array<std::int64_t, 3> long_address{};
  for (int j = 0; j < 3; ++j) {
    long_address[static_cast<std::size_t>(j)] =
        static_cast<std::int64_t>(address_double[j]) *
        divisor[static_cast<std::size_t>(j)];
  }

  std::size_t best = self;
  for (auto const &rot : rot_reciprocal) {
    Vector3i rotated;
    bool indivisible = false;
    for (int k = 0; k < 3; ++k) {
      std::int64_t const value =
          static_cast<std::int64_t>(rot(k, 0)) * long_address[0] +
          static_cast<std::int64_t>(rot(k, 1)) * long_address[1] +
          static_cast<std::int64_t>(rot(k, 2)) * long_address[2];
      auto const div = divisor[static_cast<std::size_t>(k)];
      if (value % div != 0) {
        indivisible = true;
        break;
      }
      rotated[k] = static_cast<int>(value / div);
      bool const odd = rotated[k] % 2 != 0;
      if ((odd && is_shift[k] == 0) || (!odd && is_shift[k] == 1)) {
        indivisible = true;
        break;
      }
    }
    if (!indivisible) {
      best = std::min(best, grid_point_double_mesh(rotated, mesh));
    }
  }
  return best;
}

} // namespace

IrReciprocalMesh
ir_reciprocal_mesh(Vector3i const &mesh, Vector3i const &is_shift,
                   std::vector<Matrix3i> const &rot_reciprocal) {
  std::vector<Vector3i> grid_address = all_grid_addresses(mesh);
  std::size_t const n = grid_address.size();
  bool const normal =
      mesh_has_conventional_symmetry(mesh, is_shift, rot_reciprocal);

  std::array<std::int64_t, 3> divisor{};
  if (!normal) {
    divisor = std::array<std::int64_t, 3>{
        static_cast<std::int64_t>(mesh[1]) * mesh[2],
        static_cast<std::int64_t>(mesh[2]) * mesh[0],
        static_cast<std::int64_t>(mesh[0]) * mesh[1],
    };
  }

  std::vector<std::size_t> mapping(n);
  std::ranges::transform(
      grid_address | std::views::enumerate, mapping.begin(),
      [&](auto const &entry) {
        auto const &[index, address] = entry;
        auto const i = static_cast<std::size_t>(index);
        Vector3i const address_double =
            address_double_mesh(address, mesh, is_shift);
        return normal ? normal_representative(address_double, mesh, i,
                                              rot_reciprocal)
                      : distortion_representative(address_double, mesh,
                                                  is_shift, i, divisor,
                                                  rot_reciprocal);
      });

  auto const num_ir = static_cast<std::size_t>(std::ranges::count_if(
      std::views::iota(std::size_t{0}, n),
      [&](std::size_t i) { return mapping[i] == i; }));
  return {std::move(grid_address), std::move(mapping), num_ir};
}

Result<IrReciprocalMesh>
ir_reciprocal_mesh(Cell const &cell, Vector3i const &mesh,
                   Vector3i const &is_shift, TimeReversal time_reversal,
                   Tolerance const &tol) {
  if (mesh[0] <= 0 || mesh[1] <= 0 || mesh[2] <= 0) {
    return leaf::new_error(e_invalid_mesh{});
  }
  BOOST_LEAF_AUTO(dataset,
                  analysis::SymmetryAnalyzer::from_cell(cell, tol).dataset());
  std::vector<Matrix3i> rotations;
  rotations.reserve(dataset.operations.size());
  std::ranges::transform(dataset.operations, std::back_inserter(rotations),
                         [](auto const &op) { return op.rotation; });
  return ir_reciprocal_mesh(mesh, is_shift,
                            point_group_reciprocal(rotations, time_reversal));
}

IrReciprocalMesh
stabilized_reciprocal_mesh(Vector3i const &mesh, Vector3i const &is_shift,
                           TimeReversal time_reversal,
                           std::vector<Matrix3i> const &rotations,
                           std::vector<Vector3d> const &qpoints) {
  auto const rot_reciprocal = point_group_reciprocal(rotations, time_reversal);
  const double tolerance = 0.01 / mesh.sum();
  auto const rot_q =
      point_group_reciprocal_with_q(rot_reciprocal, tolerance, qpoints);
  return ir_reciprocal_mesh(mesh, is_shift, rot_q);
}

std::vector<std::size_t>
grid_points_by_rotations(Vector3i const &address_orig,
                         std::vector<Matrix3i> const &rot_reciprocal,
                         Vector3i const &mesh, Vector3i const &is_shift) {
  Vector3i address_double_orig = 2 * address_orig.array() + is_shift.array();
  std::vector<std::size_t> out;
  out.reserve(rot_reciprocal.size());
  std::ranges::transform(
      rot_reciprocal, std::back_inserter(out), [&](Matrix3i const &rot) {
        return grid_point_double_mesh(Vector3i(rot * address_double_orig), mesh);
      });
  return out;
}

} // namespace cppcrystal::kpoint
