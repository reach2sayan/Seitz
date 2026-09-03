#pragma once

#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/types.hpp>

namespace cppcrystal {

namespace data {
struct SpacegroupType;
}

// The matched space group of a cell or an operation set: the Hall setting it
// was matched to, the conventional/bravais lattice it was matched in, and the
// origin shift that aligns the operations with the Hall-symbol database. The
// metadata itself is not copied — it is one lookup away through the key, which
// is what type() does. (Out of line: the catalog that answers it is built on
// top of this header.)
struct SpacegroupMatch {
  HallNumber hall;
  Matrix3d bravais_lattice{Matrix3d::Identity()};
  Vector3d origin_shift{Vector3d::Zero()};

  [[nodiscard]] data::SpacegroupType const &type() const noexcept;
};

// Whether the lattice handed to OperationSet::spacegroup is the conventional
// cell (recover the primitive setting via the transformation the operations
// imply) or already a primitive cell. A template argument, so the branch is
// resolved at compile time.
enum class LatticeSetting { conventional, primitive };

} // namespace cppcrystal
