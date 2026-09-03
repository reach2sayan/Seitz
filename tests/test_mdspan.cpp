#include <cppcrystal/core/mdspan.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <vector>

using namespace cppcrystal;

namespace {
// A 2 x 3 table over a flat constexpr array is itself a constant expression.
constexpr std::array<int, 6> kFlat{0, 1, 2, 3, 4, 5};
constexpr md::table<int, 2, 3> kTable(kFlat.data());
static_assert(kTable.extent(0) == 2 && kTable.extent(1) == 3);
static_assert(kTable[1, 2] == 5);
static_assert(kTable[0, 1] == 1);
} // namespace

TEST_CASE("matrix_view is a row-major view over a flat buffer", "[mdspan]") {
  std::vector<int> buffer{10, 11, 12, 20, 21, 22};
  md::matrix_view<int> rows(buffer.data(), 2, 3);
  CHECK(rows[0, 0] == 10);
  CHECK(rows[1, 2] == 22);
  rows[1, 0] = 99;
  CHECK(buffer[3] == 99);

  // A row slice is a 1-D view sharing the buffer.
  auto const row1 = md::submdspan(rows, 1, md::full_extent);
  CHECK(row1.extent(0) == 3);
  CHECK(row1[0] == 99);
}
