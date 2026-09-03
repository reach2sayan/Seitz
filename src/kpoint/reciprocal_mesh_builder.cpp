#include <cppcrystal/kpoint/reciprocal_mesh_builder.hpp>

#include <utility>

namespace cppcrystal::kpoint {

ReciprocalMeshBuilder ReciprocalMeshBuilder::from_cell(Cell cell,
                                                       Tolerance tol) {
  return ReciprocalMeshBuilder{std::move(cell), tol};
}

Result<IrReciprocalMesh>
ReciprocalMeshBuilder::irreducible(Vector3i const &mesh,
                                   Vector3i const &is_shift,
                                   TimeReversal time_reversal) const {
  return ir_reciprocal_mesh(cell_, mesh, is_shift, time_reversal, tol_);
}

} // namespace cppcrystal::kpoint
