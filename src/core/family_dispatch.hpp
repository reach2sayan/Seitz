#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/periodicity.hpp>

#include <concepts>
#include <utility>

namespace seitz {

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

// The ONE runtime branch on the group family. Everything below it is templated
// on GroupFamily, so the family-dependent steps of the pipeline are resolved at
// compile time; this is where a cell's periodicity is turned into that
// parameter, once, at the top.
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
