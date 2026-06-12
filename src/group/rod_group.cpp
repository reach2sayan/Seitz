#include <spglib/group/rod_group.hpp>

#include <spglib/data/rod_database.hpp>

#include <span>
#include <vector>

namespace spglib::group {

namespace {

// The identity operation of an operation set (the rod group always contains it).
[[nodiscard]] SymmetryOperation identity_of(SymmetryOperations const &ops) {
  for (auto const &op : ops) {
    if (op.is_identity_rotation() &&
        op.translation.cwiseAbs().maxCoeff() < 1e-9) {
      return op;
    }
  }
  return SymmetryOperation{}; // fallback: default identity
}

} // namespace

Vector3d RodWyckoff::sample(std::span<double const> params) const {
  Vector3d p = locus_origin_;
  for (int i = 0; i < dof_ && i < static_cast<int>(params.size()); ++i) {
    p += params[static_cast<std::size_t>(i)] * locus_basis_.col(i);
  }
  return p;
}

Result<RodGroup> RodGroup::from_number(int number) {
  if (number < 1 || number > data::num_rod_groups()) {
    return leaf::new_error(e_message{
        "RodGroup::from_number: rod-group number out of range "
        "(expected 1..75)"});
  }

  SymmetryOperations ops = data::rod_operations_from_database(number);
  if (ops.empty()) {
    return leaf::new_error(e_message{
        "RodGroup::from_number: no operations for this rod group (is the "
        "generated rod table present? run tools/transcribe_rod_groups.py)"});
  }

  RodGroup rg;
  rg.number_ = number;
  rg.symbol_ = data::rod_symbol(number);
  rg.periodic_axis_ = data::rod_periodic_axis();
  rg.operations_ = std::move(ops);

  // --- General position (always correct, no table beyond the operations) ------
  // A generic point's stabiliser is trivial, so every operation is its own coset
  // representative and the multiplicity equals the group order. The locus is the
  // whole cell (origin 0, the three Cartesian basis directions), dof 3.
  {
    SymmetryOperations const orbit_ops = rg.operations_;
    SymmetryOperations const site_symmetry{identity_of(rg.operations_)};
    rg.positions_.push_back(RodWyckoff{static_cast<int>(rg.operations_.size()),
                                       /*dof=*/3, /*letter=*/'a',
                                       Vector3d::Zero(), Matrix3d::Identity(),
                                       orbit_ops, site_symmetry});
  }

  // --- Special positions (TODO) -----------------------------------------------
  // The remaining Wyckoff positions are the orbit types of the affine fixed-locus
  // arrangement: for each operation h(p) = R p + t, its fixed locus solves
  // (R - I) p == -t (exactly along the two aperiodic axes; modulo 1 along the
  // periodic axis), an affine subspace `offset + ker(R - I)`. Screw/glide
  // operations whose fixed equation is inconsistent (along an aperiodic axis)
  // contribute no locus. Close the loci under intersection, group them into
  // orbits under the group, and per orbit emit a RodWyckoff exactly as
  // group::PointGroup does for the linear (origin-fixed) 0D case — generalised
  // here to the affine, translation-bearing setting. This needs the real rod
  // table to validate (op-invariance round-trip), so it is deferred until
  // data/rod_group_tables.hpp is generated. Until then only the general position
  // is available, which is sufficient for general_position_only generation.
  //
  // When added, collect all positions, sort by ascending multiplicity (then
  // dof), and reassign letters 'a'.. (the general position ending up last), the
  // same ordering group::PointGroup uses.

  return rg;
}

} // namespace spglib::group
