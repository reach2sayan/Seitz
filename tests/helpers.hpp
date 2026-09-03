#pragma once

#include <cppcrystal/core/error.hpp>
#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/core/keys.hpp>

#include <catch2/catch_test_macros.hpp>

#include <stdexcept>
#include <utility>

// Shared Result<T> plumbing for the test suites, so each file does not carry
// its own copy of the same two try_handle_all wrappers.
namespace cppcrystal::test {

// Shorthand for the validated keys: the tests name literal settings by the
// dozen, and `*HallNumber::of(GroupFamily::space, 1)` at each of them buries
// what is being tested.
[[nodiscard]] constexpr HallNumber space_hall(int index) noexcept {
  return *HallNumber::of(GroupFamily::space, index);
}
[[nodiscard]] constexpr HallNumber layer_hall(int index) noexcept {
  return *HallNumber::of(GroupFamily::layer, index);
}
[[nodiscard]] constexpr UniNumber uni_number(int number) noexcept {
  return *UniNumber::of(number);
}

// The determination of one cell, the way callers reach it now. A shorthand
// only: SymmetryAnalyzer::from_cell(...).dataset() at two dozen call sites
// would bury what each test is actually checking.
[[nodiscard]] inline Result<analysis::Dataset>
dataset_of(Cell const &cell, Tolerance const &tol = {}) {
  return analysis::SymmetryAnalyzer::from_cell(cell, tol).dataset();
}

[[nodiscard]] inline Result<analysis::MagneticDataset>
magnetic_dataset_of(MagneticCell const &cell, MagneticTolerance const &tol = {}) {
  return analysis::MagneticSymmetryAnalyzer::from_cell(cell, tol).dataset();
}

// The success value of `r`, failing the current test if it carries an error.
template <class T> T must(Result<T> r) {
  return leaf::try_handle_all(
      [&]() -> Result<T> { return std::move(r); },
      [](leaf::error_info const &) -> T {
        FAIL("unexpected error result");
        throw std::logic_error("unreachable");
      });
}

// Whether `make()` produced an error — the expected-failure counterpart of
// must(). Takes the call rather than the Result so a BOOST_LEAF_AUTO chain can
// be written inline.
template <class F> [[nodiscard]] bool errored(F &&make) {
  return leaf::try_handle_all(
      [&]() -> Result<bool> {
        if (auto r = std::forward<F>(make)(); !r) {
          return r.error();
        }
        return false;
      },
      [](leaf::error_info const &) { return true; });
}

} // namespace cppcrystal::test
