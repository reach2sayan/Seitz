#include <cppcrystal/analysis/operation_set.hpp>

namespace cppcrystal::analysis {

std::vector<Matrix3i> OperationSet::rotations() const {
  return symmetry::rotations_of(ops_);
}

std::optional<std::pair<SymmetryOperations, Matrix3d>>
OperationSet::to_primitive(double symprec) const {
  return symmetry::primitive_symmetry(ops_, symprec);
}

Result<spacegroup::Spacegroup>
OperationSet::spacegroup_type(Matrix3d const &lattice, double symprec,
                              spacegroup::LatticeSetting setting) const {
  using spacegroup::LatticeSetting;
  if (setting == LatticeSetting::primitive) {
    return spacegroup::spacegroup_type_from_symmetry<LatticeSetting::primitive>(
        ops_, lattice, symprec);
  }
  return spacegroup::spacegroup_type_from_symmetry<LatticeSetting::conventional>(
      ops_, lattice, symprec);
}

Result<spacegroup::Spacegroup>
OperationSet::search_spacegroup(Matrix3d const &prim_lattice,
                                double symprec) const {
  return spacegroup::search_spacegroup_with_symmetry(ops_, prim_lattice,
                                                     symprec);
}

std::vector<Vector3d> MagneticOperationSet::pure_translations() const {
  return spin::collect_pure_translations(ops_);
}

Result<magnetic::MagneticTypeIdentification>
MagneticOperationSet::identify_spacegroup_type(Matrix3d const &lattice,
                                               double symprec) const {
  return magnetic::identify_magnetic_spacegroup_type(lattice, ops_, symprec);
}

} // namespace cppcrystal::analysis
