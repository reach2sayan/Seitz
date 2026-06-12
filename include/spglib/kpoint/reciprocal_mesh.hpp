#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>

#include <cstddef>
#include <vector>

// Irreducible reciprocal-space sampling (the symmetry reduction of a k-point
// mesh). Port of spglib's kpoint.c (3D path). Reciprocal-space rotations are the
// transpose of the real-space rotations. This header grows over Phase 8: the
// reciprocal point group lands first; the irreducible / stabilized meshes follow.
namespace spglib::kpoint {

// The reciprocal-space point group implied by a set of real-space rotations:
// each reciprocal rotation is the transpose of a real-space rotation; with time
// reversal the inversion partner (-transpose) is added as well. Duplicates are
// removed, keeping first occurrence. Port of kpoint.c get_point_group_reciprocal.
[[nodiscard]] std::vector<Matrix3i>
point_group_reciprocal(std::vector<Matrix3i> const &rotations,
                       bool time_reversal);

// The subset of `rot_reciprocal` that maps the q-point set onto itself (modulo a
// reciprocal lattice vector) within `symprec` — the stabilizer of the q-points.
// Port of kpoint.c get_point_group_reciprocal_with_q.
[[nodiscard]] std::vector<Matrix3i>
point_group_reciprocal_with_q(std::vector<Matrix3i> const &rot_reciprocal,
                              double symprec,
                              std::vector<Vector3d> const &qpoints);

// The irreducible reciprocal mesh: the reducible grid points and, for each, the
// grid-point index of its irreducible representative (the smallest index in its
// orbit under the reciprocal point group). `ir_mapping_table[i] == i` ⟺ i is a
// representative; `num_ir` is the number of representatives.
struct IrReciprocalMesh {
  std::vector<Vector3i> grid_address;
  std::vector<std::size_t> ir_mapping_table;
  std::size_t num_ir = 0;
};

// Irreducible mesh from a reciprocal point group (the rotation-set form, port of
// kpt_get_dense_irreducible_reciprocal_mesh). `rot_reciprocal` must already be
// the reciprocal group (time reversal, if wanted, baked in by the caller).
[[nodiscard]] IrReciprocalMesh
ir_reciprocal_mesh(Vector3i const &mesh, Vector3i const &is_shift,
                   std::vector<Matrix3i> const &rot_reciprocal);

// Irreducible mesh of a crystal: derive the space-group rotations via
// get_dataset, build the reciprocal point group (with time reversal if
// requested), then reduce the mesh. Port of spglib.c get_(dense_)ir_reciprocal_
// mesh. Errors via get_dataset (e_spacegroup_search_failed / …).
// For the object-oriented entry point, prefer
// spglib::kpoint::ReciprocalMeshBuilder::irreducible(...).
[[nodiscard]] Result<IrReciprocalMesh>
ir_reciprocal_mesh(Cell const &cell, Vector3i const &mesh,
                   Vector3i const &is_shift, bool time_reversal,
                   double symprec = kDefaultSymprec,
                   AngleTolerance angle_tolerance = std::nullopt);

// Irreducible mesh stabilized by a set of q-points: reduce only by the
// rotations that map the q-point set onto itself. Port of
// kpt_get_dense_stabilized_reciprocal_mesh.
[[nodiscard]] IrReciprocalMesh
stabilized_reciprocal_mesh(Vector3i const &mesh, Vector3i const &is_shift,
                           bool time_reversal,
                           std::vector<Matrix3i> const &rotations,
                           std::vector<Vector3d> const &qpoints);

// The grid-point indices that one grid address maps to under each reciprocal
// rotation. Port of kpt_get_dense_grid_points_by_rotations.
[[nodiscard]] std::vector<std::size_t>
grid_points_by_rotations(Vector3i const &address_orig,
                         std::vector<Matrix3i> const &rot_reciprocal,
                         Vector3i const &mesh, Vector3i const &is_shift);

} // namespace spglib::kpoint
