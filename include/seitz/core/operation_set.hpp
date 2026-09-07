#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/fractional.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/lattice.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/spacegroup_match.hpp>

#include <boost/container/small_vector.hpp>

#include <concepts>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz {

// One element's spatial part, and whether the family carries time reversal:
// what lets OperationSet be written once for both families.
template <class Op> struct OperationTraits;

template <> struct OperationTraits<SymmetryOperation> {
  static constexpr bool has_time_reversal = false;
  [[nodiscard]] static constexpr SymmetryOperation const &
  spatial(SymmetryOperation const &op) noexcept {
    return op;
  }
};

template <> struct OperationTraits<MagneticSymmetryOperation> {
  static constexpr bool has_time_reversal = true;
  [[nodiscard]] static constexpr SymmetryOperation const &
  spatial(MagneticSymmetryOperation const &op) noexcept {
    return op.spatial;
  }
};

template <class Op>
concept Operation = requires(Op const &op) {
  {
    OperationTraits<Op>::spatial(op)
  } -> std::convertible_to<SymmetryOperation>;
  { OperationTraits<Op>::has_time_reversal } -> std::convertible_to<bool>;
};

template <Operation Op> class OperationSet;

namespace detail {

// The primitive operations and the primitive->conventional transformation a
// conventional set implies. Defined in src/core, so to_primitive stays inline
// without this header naming a pipeline type.
[[nodiscard]] std::optional<std::pair<std::vector<SymmetryOperation>, Matrix3d>>
primitive_operations(std::span<SymmetryOperation const> operations,
                     Tolerance const &tol);

// The space group these operations imply in `lattice`; the matcher is private
// to src/spacegroup.
[[nodiscard]] Result<SpacegroupMatch>
spacegroup_of_operations(OperationSet<SymmetryOperation> const &operations,
                         Lattice const &lattice, LatticeSetting setting,
                         Tolerance const &tol);

// The point group of a rotation set; the identifier is private to
// src/symmetry. `layer_axis` (a layer cell's aperiodic axis in the input
// basis) is data the layer family sorts to c and the space family ignores.
template <GroupFamily F>
[[nodiscard]] Result<PointGroupMatch>
point_group_of_rotations(std::span<Matrix3i const> rotations,
                         std::optional<int> layer_axis);
extern template Result<PointGroupMatch>
    point_group_of_rotations<GroupFamily::space>(std::span<Matrix3i const>,
                                                 std::optional<int>);
extern template Result<PointGroupMatch>
    point_group_of_rotations<GroupFamily::layer>(std::span<Matrix3i const>,
                                                 std::optional<int>);

// The magnetic space group these operations imply in `lattice`; the
// identifier is private to src/magnetic.
[[nodiscard]] Result<MagneticMatch> magnetic_spacegroup_of_operations(
    OperationSet<MagneticSymmetryOperation> const &operations,
    Lattice const &lattice, Tolerance const &tol);

} // namespace detail

// An immutable set of symmetry operations. A range, so it composes with
// std::views; built once and never mutated.
template <Operation Op> class OperationSet {
public:
  using value_type = Op;
  using Traits = OperationTraits<Op>;

  OperationSet() = default;
  explicit OperationSet(std::vector<Op> ops) noexcept : ops_(std::move(ops)) {}

  // From any range, so a filtered/transformed view materialises into a set.
  template <std::ranges::input_range R>
    requires std::constructible_from<Op, std::ranges::range_reference_t<R>>
  explicit OperationSet(std::from_range_t, R &&range)
      : ops_(std::ranges::begin(range), std::ranges::end(range)) {}

  [[nodiscard]] auto begin() const noexcept { return ops_.begin(); }
  [[nodiscard]] auto end() const noexcept { return ops_.end(); }
  [[nodiscard]] std::size_t size() const noexcept { return ops_.size(); }
  [[nodiscard]] bool empty() const noexcept { return ops_.empty(); }
  [[nodiscard]] Op const &operator[](std::size_t i) const noexcept {
    return ops_[i];
  }
  [[nodiscard]] std::span<Op const> span() const noexcept { return ops_; }

  // The rotation parts in order; de-duplication is the caller's business.
  [[nodiscard]] std::vector<Matrix3i> rotations() const {
    return {std::from_range, ops_ | std::views::transform([](Op const &op) {
                               return Traits::spatial(op).rotation;
                             })};
  }

  // Translations of the identity-rotation ops, anti-translations excluded;
  // the zero translation included.
  [[nodiscard]] std::vector<Vector3d> pure_translations() const {
    // Not const: filter_view caches its begin(), so it is not const-iterable.
    auto pure = ops_ | std::views::filter([](Op const &op) {
                  if constexpr (Traits::has_time_reversal) {
                    if (op.time_reversal) {
                      return false;
                    }
                  }
                  return Traits::spatial(op).is_identity_rotation();
                }) |
                std::views::transform([](Op const &op) {
                  return Traits::spatial(op).translation;
                });
    return {std::from_range, pure};
  }

  // (T, 0)(R, t)(T, 0)^-1 for every op; a time-reversal flag rides along.
  [[nodiscard]] OperationSet conjugated_by(Matrix3d const &t,
                                           Matrix3d const &t_inv) const {
    // The formula lives once, in symmetry_operation.hpp; applied to the
    // spatial part only.
    return OperationSet{std::from_range,
                        ops_ | std::views::transform([&](Op op) {
                          SymmetryOperation &sp = spatial_of(op);
                          sp = seitz::conjugated_by(sp, t, t_inv);
                          return op;
                        })};
  }

  // The underlying space-group operations, dropping the time-reversal flags.
  [[nodiscard]] OperationSet<SymmetryOperation> spatial() const
    requires Traits::has_time_reversal
  {
    return OperationSet<SymmetryOperation>{
        std::from_range, ops_ | std::views::transform([](Op const &op) {
                           return Traits::spatial(op);
                         })};
  }

  // The space group these operations imply, with no atomic positions: build
  // the primitive symmetry, Niggli-reduce the implied primitive lattice, match
  // it against the Hall database. `S` says whether `lattice` is conventional or
  // already primitive. Errors with e_spacegroup_search_failed.
  template <LatticeSetting S = LatticeSetting::conventional>
  [[nodiscard]] Result<SpacegroupMatch> spacegroup(Lattice const &lattice,
                                                   Tolerance const &tol) const
    requires(!Traits::has_time_reversal)
  {
    return detail::spacegroup_of_operations(*this, lattice, S, tol);
  }

  // The magnetic space group these operations imply in `lattice`: the UNI
  // number, its type, the reference space group and the standardizing
  // setting. Errors with e_magnetic_symmetry_search_failed.
  [[nodiscard]] Result<MagneticMatch> spacegroup(Lattice const &lattice,
                                                 Tolerance const &tol) const
    requires Traits::has_time_reversal
  {
    return detail::magnetic_spacegroup_of_operations(*this, lattice, tol);
  }

  // The point group of the rotation parts, with the change of basis to its
  // conventional axes. The layer family rejects the cubic groups and sorts
  // `layer_axis` to c. Errors with e_pointgroup_not_found.
  template <GroupFamily F = GroupFamily::space>
  [[nodiscard]] Result<PointGroupMatch>
  point_group(std::optional<int> layer_axis = std::nullopt) const {
    // A primitive set has at most 48 rotations: no heap unless centred.
    auto const view = ops_ | std::views::transform([](Op const &op) {
                        return Traits::spatial(op).rotation;
                      });
    boost::container::small_vector<Matrix3i, 48> const rotations(
        std::ranges::begin(view), std::ranges::end(view));
    return detail::point_group_of_rotations<F>(
        std::span{rotations.data(), rotations.size()}, layer_axis);
  }

  // The primitive operations these (conventional) ones imply, with the
  // primitive -> conventional transformation t_mat:
  // (a_p, b_p, c_p) . t_mat = (a_c, b_c, c_c). The pure translations give the
  // primitive lattice in translation space; the distinct rotations, carried to
  // that setting, give the operations. nullopt if the set is inconsistent.
  [[nodiscard]] std::optional<std::pair<OperationSet, Matrix3d>>
  to_primitive(Tolerance const &tol) const
    requires(!Traits::has_time_reversal)
  {
    auto primitive = detail::primitive_operations(span(), tol);
    if (!primitive) {
      return std::nullopt;
    }
    return std::pair{OperationSet{std::move(primitive->first)},
                     primitive->second};
  }

private:
  [[nodiscard]] static SymmetryOperation &spatial_of(Op &op) noexcept {
    if constexpr (Traits::has_time_reversal) {
      return op.spatial;
    } else {
      return op;
    }
  }

  std::vector<Op> ops_;
};

using Operations = OperationSet<SymmetryOperation>;
using MagneticOperations = OperationSet<MagneticSymmetryOperation>;

} // namespace seitz

#pragma GCC visibility pop
