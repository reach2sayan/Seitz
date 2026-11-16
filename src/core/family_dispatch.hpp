#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/periodicity.hpp>

#include <concepts>
#include <utility>

namespace seitz {

// In the conventional/standardized setting a layer group always has its
// aperiodic axis as c (axis 2); a 3D group is periodic on all three.
template <GroupFamily F>
[[nodiscard]] constexpr CellPeriodicity conventional_periodicity() noexcept {
  if constexpr (F == GroupFamily::layer) {
    return aperiodic_along(2);
  } else {
    return all_periodic();
  }
}

template <class Fn, GroupFamily F>
using family_call_result_t =
    decltype(std::declval<Fn>().template operator()<F>());

// A visitor for dispatch_family: callable with the family as an explicit
// template argument, and agreeing on the return type across both families
// (dispatch_family returns a ternary over the two, so a disagreement has to
// be diagnosed here rather than inside the body).
template <class Fn>
concept FamilyVisitor =
    requires(Fn &&fn) {
      std::forward<Fn>(fn).template operator()<GroupFamily::layer>();
      std::forward<Fn>(fn).template operator()<GroupFamily::space>();
    } && std::same_as<family_call_result_t<Fn, GroupFamily::layer>,
                      family_call_result_t<Fn, GroupFamily::space>>;

// The ONE runtime branch on the group family: everything below is templated on
// GroupFamily, so the family-dependent steps resolve at compile time. Here a
// cell's periodicity becomes that parameter, once, at the top.
//
//   dispatch_family(cell.periodicity(),
//                   [&]<GroupFamily F>() { return run_pipeline<F>(); });
template <FamilyVisitor Fn>
[[nodiscard]] decltype(auto) dispatch_family(CellPeriodicity const &periodicity,
                                             Fn &&fn) {
  return family_of(periodicity) == GroupFamily::layer
             ? std::forward<Fn>(fn).template operator()<GroupFamily::layer>()
             : std::forward<Fn>(fn).template operator()<GroupFamily::space>();
}

} // namespace seitz
