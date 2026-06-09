#include <spglib/core/error.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace spglib;
using Catch::Approx;

namespace {
Result<int> may_fail(bool ok) {
  if (!ok)
    return leaf::new_error(e_atoms_too_close{0.3});
  return 42;
}
} // namespace

TEST_CASE("Result carries a value on success", "[error]") {
  auto r = may_fail(true);
  REQUIRE(r);
  CHECK(r.value() == 42);
}

TEST_CASE("Result propagates an error tag through BOOST_LEAF_CHECK",
          "[error]") {
  double seen = -1.0;
  int code = leaf::try_handle_all(
      [&]() -> Result<int> {
        int v = BOOST_LEAF_CHECK(may_fail(false));
        return v;
      },
      [&](e_atoms_too_close const &e) {
        seen = e.distance;
        return 1;
      },
      [] { return 2; });
  CHECK(code == 1);
  CHECK(seen == Approx(0.3));
}
