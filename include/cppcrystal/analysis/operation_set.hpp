#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_symmetry_operation.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>
#include <cppcrystal/magnetic/magnetic_spacegroup.hpp>
#include <cppcrystal/spacegroup/spacegroup.hpp>
#include <cppcrystal/spin/spin.hpp>
#include <cppcrystal/symmetry/pointgroup.hpp>
#include <cppcrystal/symmetry/primitive.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace cppcrystal::analysis {

// CRTP base of the operation-set wrappers: the shared construction and
// container face over an owned operations vector. Static polymorphism only —
// `Derived` is a concrete type, there is no runtime hierarchy.
template <class Derived, class Ops> class OperationSetBase {
public:
  [[nodiscard]] static Derived from_operations(Ops ops) {
    return Derived{std::move(ops)};
  }
  [[nodiscard]] Ops const &operations() const noexcept { return ops_; }
  [[nodiscard]] std::size_t size() const noexcept { return ops_.size(); }

protected:
  explicit OperationSetBase(Ops ops) : ops_(std::move(ops)) {}
  Ops ops_;
};

// An immutable wrapper that gives a bare SymmetryOperations vector an object
// identity: the free functions that take SymmetryOperations as their de-facto
// `this` (rotations_of, primitive_symmetry, spacegroup_type_from_symmetry,
// search_spacegroup_with_symmetry) become methods here. Nothing is memoized
// (projections are cheap, an OperationSet is usually transient).
class OperationSet : public OperationSetBase<OperationSet, SymmetryOperations> {
public:
  // The rotation parts of the operations, de-duplicated by value internally
  // (symmetry::rotations_of).
  [[nodiscard]] std::vector<Matrix3i> rotations() const;

  // The primitive operations together with the primitive->conventional
  // transformation matrix (symmetry::primitive_symmetry). nullopt if the set is
  // inconsistent (the implied primitive lattice cannot be resolved).
  [[nodiscard]] std::optional<std::pair<SymmetryOperations, Matrix3d>>
  to_primitive(double symprec) const;

  // The space group implied by these operations and a `lattice`
  // (spacegroup::spacegroup_type_from_symmetry). `setting` selects whether
  // `lattice` is the conventional cell (recover the primitive via the implied
  // t_mat) or already primitive.
  [[nodiscard]] Result<spacegroup::Spacegroup>
  spacegroup_type(Matrix3d const &lattice, double symprec,
                  spacegroup::LatticeSetting setting =
                      spacegroup::LatticeSetting::conventional) const;

  // The space group of a primitive cell given by `prim_lattice` and these
  // operations (spacegroup::search_spacegroup_with_symmetry).
  [[nodiscard]] Result<spacegroup::Spacegroup>
  search_spacegroup(Matrix3d const &prim_lattice, double symprec) const;

private:
  friend OperationSetBase;
  using OperationSetBase::OperationSetBase;
};

// The magnetic counterpart of OperationSet: an immutable wrapper over a
// MagneticSymmetryOperations vector, exposing the free functions that take it as
// their `this`.
class MagneticOperationSet
    : public OperationSetBase<MagneticOperationSet,
                              MagneticSymmetryOperations> {
public:
  // The pure (non-time-reversal) translations, including the zero translation
  // (spin::collect_pure_translations).
  [[nodiscard]] std::vector<Vector3d> pure_translations() const;

  // The magnetic space-group identification of these operations in `lattice`
  // (magnetic::identify_magnetic_spacegroup_type).
  [[nodiscard]] Result<magnetic::MagneticTypeIdentification>
  identify_spacegroup_type(Matrix3d const &lattice, double symprec) const;

private:
  friend OperationSetBase;
  using OperationSetBase::OperationSetBase;
};

} // namespace cppcrystal::analysis
