#pragma once

#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/core/periodicity.hpp>

#include <utility>

namespace cppcrystal {

// The ONE runtime branch on the group family. Everything below it is templated
// on GroupFamily, so the family-dependent steps of the pipeline are resolved at
// compile time; this is where a cell's periodicity is turned into that
// parameter, once, at the top.
//
//   dispatch_family(cell.periodicity(),
//                   [&]<GroupFamily F>() { return run_pipeline<F>(); });
template <class Fn>
[[nodiscard]] decltype(auto) dispatch_family(CellPeriodicity const &periodicity,
                                             Fn &&fn) {
  return family_of(periodicity) == GroupFamily::layer
             ? std::forward<Fn>(fn).template operator()<GroupFamily::layer>()
             : std::forward<Fn>(fn).template operator()<GroupFamily::space>();
}

} // namespace cppcrystal
