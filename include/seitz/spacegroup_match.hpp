#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/point_group.hpp>
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

// The matched point group of a rotation set, with the integer change of basis
// that brings the rotations into the conventional setting (columns are the
// chosen axes).
struct PointGroupMatch {
  PointGroup type;
  Matrix3i transformation{Matrix3i::Zero()};
};

// input cell <--> standardized setting:
// (a) the change of basis + origin shift to align ops with the database,
// (b) rigid rotation to the idealized standardized lattice.
struct Setting {
  Matrix3d transformation{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};
  Matrix3d rigid_rotation{Matrix3d::Identity()};
};

// Construction type of a magnetic space group (Barnighausen / BNS types I-IV).
enum class MagneticType { type_i = 1, type_ii = 2, type_iii = 3, type_iv = 4 };

// The matched magnetic space group of a magnetic operation set: the UNI
// number and its type, the family (types I-III) or maximal (type IV) space
// group, and the setting that standardizes it.
struct MagneticMatch {
  UniNumber uni;
  MagneticType type = MagneticType::type_i;
  HallNumber hall;
  Setting setting;
};

} // namespace seitz

#pragma GCC visibility pop
