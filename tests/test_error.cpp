#include <seitz/core/error.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace seitz;
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
        // BOOST_LEAF_AUTO, not `v = BOOST_LEAF_CHECK(...)`: LEAF only expands
        // BOOST_LEAF_CHECK to a value-yielding expression where GNU statement
        // expressions exist. Under MSVC it is a statement, so the value form
        // does not compile. BOOST_LEAF_AUTO is the portable spelling.
        BOOST_LEAF_AUTO(v, may_fail(false));
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
