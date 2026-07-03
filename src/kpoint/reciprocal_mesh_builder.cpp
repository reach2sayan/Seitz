#include <cppcrystal/kpoint/reciprocal_mesh_builder.hpp>

#include <utility>

namespace cppcrystal::kpoint {

ReciprocalMeshBuilder ReciprocalMeshBuilder::from_cell(
    Cell cell, double symprec, AngleTolerance angle_tolerance) {
  return ReciprocalMeshBuilder{std::move(cell),
                               Tolerance{symprec, angle_tolerance}};
}

Result<IrReciprocalMesh> ReciprocalMeshBuilder::irreducible(
    Vector3i const &mesh, Vector3i const &is_shift, bool time_reversal) const {
  return ir_reciprocal_mesh(cell_, mesh, is_shift, time_reversal, tol_.symprec,
                            tol_.angle_tolerance);
}

} // namespace cppcrystal::kpoint
