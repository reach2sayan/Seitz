#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/types.hpp>

#pragma GCC visibility push(default)

namespace seitz {

namespace data {
struct SpacegroupType;
}

// The matched space group of a cell or operation set: the Hall setting, the
// conventional (bravais) lattice it was matched in, and the origin shift
// aligning the operations with the Hall-symbol database. The metadata is one
// key lookup away -- type(), out of line, since its catalog is built on this
// header.
struct SpacegroupMatch {
  HallNumber hall;
  Matrix3d bravais_lattice{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};
  [[nodiscard]] data::SpacegroupType const &type() const noexcept;
};

// Whether the lattice given to OperationSet::spacegroup is conventional (the
// primitive setting is recovered from the transformation the operations imply)
// or already primitive.
enum class LatticeSetting { conventional, primitive };

} // namespace seitz

#pragma GCC visibility pop
