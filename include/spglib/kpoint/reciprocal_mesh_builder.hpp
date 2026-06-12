#pragma once

#include <spglib/core/cell.hpp>
#include <spglib/core/error.hpp>
#include <spglib/core/tolerance.hpp>
#include <spglib/core/types.hpp>
#include <spglib/kpoint/reciprocal_mesh.hpp>

namespace spglib::kpoint {

// Object-oriented entry point for the crystal irreducible-mesh reduction: a
// persistent view over a Cell + tolerances that builds the irreducible
// reciprocal mesh. The functional core (ir_reciprocal_mesh) is left untouched;
// this facade owns the inputs once so a caller can reduce several meshes of the
// same crystal without re-threading the cell and tolerances. Lives in kpoint/
// (not on SymmetryAnalyzer) so the IrReciprocalMesh / kpoint includes stay out
// of the core analysis header.
//
// Immutable after construction; no setters. Unlike SymmetryAnalyzer this does
// not memoize — each call takes a different mesh, so there is nothing to cache.
class ReciprocalMeshBuilder {
public:
  [[nodiscard]] static ReciprocalMeshBuilder
  from_cell(Cell cell, double symprec = kDefaultSymprec,
            AngleTolerance angle_tolerance = std::nullopt);

  [[nodiscard]] Cell const &cell() const noexcept { return cell_; }
  [[nodiscard]] double symprec() const noexcept { return tol_.symprec; }
  [[nodiscard]] AngleTolerance angle_tolerance() const noexcept {
    return tol_.angle_tolerance;
  }

  // The irreducible reciprocal mesh of the crystal: derives the space-group
  // rotations via get_dataset, builds the reciprocal point group (adding the
  // inversion partner when `time_reversal`), then reduces the mesh. Forwards to
  // kpoint::ir_reciprocal_mesh(Cell, ...). Errors via get_dataset.
  [[nodiscard]] Result<IrReciprocalMesh>
  irreducible(Vector3i const &mesh, Vector3i const &is_shift,
              bool time_reversal) const;

private:
  ReciprocalMeshBuilder(Cell cell, Tolerance tol)
      : cell_(std::move(cell)), tol_(tol) {}

  Cell cell_;
  Tolerance tol_;
};

} // namespace spglib::kpoint
