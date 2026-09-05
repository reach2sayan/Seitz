#pragma once

#include <seitz/core/error.hpp>
#include <seitz/core/fractional.hpp>
#include <seitz/core/magnetic_symmetry_operation.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/tolerance.hpp>
#include <seitz/core/types.hpp>
#include <seitz/spacegroup_match.hpp>

#include <concepts>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz {

namespace detail {

// The primitive operations and the primitive->conventional transformation
// implied by a conventional operation set. Defined in src/core, over the
// private primitive-cell machinery; declared here so OperationSet::to_primitive
// can stay inline without the public header naming a pipeline type.
[[nodiscard]] std::optional<std::pair<std::vector<SymmetryOperation>, Matrix3d>>
primitive_operations(std::span<SymmetryOperation const> operations,
                     Tolerance const &tol);

// The space group these operations imply in `lattice`. Same arrangement: the
// matcher itself is private to src/spacegroup.
[[nodiscard]] Result<SpacegroupMatch>
spacegroup_of_operations(std::span<SymmetryOperation const> operations,
                         Matrix3d const &lattice, LatticeSetting setting,
                         Tolerance const &tol);

} // namespace detail

// What a set of operations can do with one element: reach its spatial part,
// and whether the family carries time reversal at all. Lets OperationSet be
// written once for both families instead of twice.
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

// An immutable set of symmetry operations: the subject that the free functions
// taking a bare operations vector used to pass around. A range, so it composes
// with std::views; built once from a vector and never mutated.
template <Operation Op> class OperationSet {
public:
  using value_type = Op;
  using Traits = OperationTraits<Op>;

  OperationSet() = default;
  explicit OperationSet(std::vector<Op> ops) noexcept : ops_(std::move(ops)) {}

  // From any range of operations, so a filtered/transformed view materialises
  // straight into a set.
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

  // The rotation parts, one per operation, in order. De-duplication is the
  // caller's business (the point-group identification does its own).
  [[nodiscard]] std::vector<Matrix3i> rotations() const {
    return {std::from_range, ops_ | std::views::transform([](Op const &op) {
                               return Traits::spatial(op).rotation;
                             })};
  }

  // The pure translations: the identity-rotation operations' translations,
  // excluding the time-reversal anti-translations. Includes the zero
  // translation.
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

  // Change of basis of every operation: (T, 0)(R, t)(T, 0)^-1. Any extra
  // fields (a time-reversal flag) ride along unchanged.
  [[nodiscard]] OperationSet conjugated_by(Matrix3d const &t,
                                           Matrix3d const &t_inv) const {
    // The per-operation formula lives once, in symmetry_operation.hpp; here it
    // is applied to the spatial part so a time-reversal flag rides along.
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

  // The primitive operations implied by these (typically conventional) ones,
  // together with the transformation `t_mat` from the primitive to the
  // conventional setting: (a_p, b_p, c_p) . t_mat = (a_c, b_c, c_c). The pure
  // translations define the primitive lattice in "translation space"; the
  // distinct rotations, transformed to that setting, are the primitive
  // operations. std::nullopt if the set is inconsistent.
  // The space group these operations imply, with no atomic positions: build
  // the primitive symmetry, Niggli-reduce the implied primitive lattice, and
  // match it against the Hall database. `S` says whether `lattice` is the
  // conventional cell or already a primitive one. Errors with
  // e_spacegroup_search_failed.
  template <LatticeSetting S = LatticeSetting::conventional>
  [[nodiscard]] Result<SpacegroupMatch> spacegroup(Matrix3d const &lattice,
                                                   Tolerance const &tol) const
    requires(!Traits::has_time_reversal)
  {
    return detail::spacegroup_of_operations(span(), lattice, S, tol);
  }

  [[nodiscard]] std::optional<std::pair<OperationSet, Matrix3d>>
  to_primitive(Tolerance const &tol) const
    requires(!Traits::has_time_reversal)
  {
    auto primitive = detail::primitive_operations(span(), tol);
    if (!primitive) {
      return std::nullopt;
    }
    return std::make_pair(OperationSet{std::move(primitive->first)},
                          primitive->second);
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
