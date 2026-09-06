#include "helpers.hpp"

#include <seitz/core/keys.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/data/spg_database.hpp>

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <ranges>
#include <string>
#include <string_view>

using namespace seitz;
using seitz::test::errored;

namespace {

// The rotation of a triplet, written the way it reads: one row per coordinate.
[[nodiscard]] Matrix3i rows(std::array<int, 9> const &r) {
  Matrix3i m;
  m << r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8];
  return m;
}

struct Reading {
  std::string_view text;
  std::array<int, 9> rotation;
  Vector3d translation;
};

} // namespace

TEST_CASE("from_xyz reads the triplet forms CIF files carry", "[xyz]") {
  auto const readings = std::array{
      Reading{"x,y,z", {1, 0, 0, 0, 1, 0, 0, 0, 1}, {0.0, 0.0, 0.0}},
      Reading{"-y,x-y,z+1/2", {0, -1, 0, 1, -1, 0, 0, 0, 1}, {0.0, 0.0, 0.5}},
      // The translation written first, which is what makes the rotation
      // reading have to backtrack over the numerator.
      Reading{"1/2+x, -y, z", {1, 0, 0, 0, -1, 0, 0, 0, 1}, {0.5, 0.0, 0.0}},
      // An explicit coefficient the grammar has to read, in a rotation that
      // is still unimodular.
      Reading{"x, 2x+y+1/2, z", {1, 0, 0, 2, 1, 0, 0, 0, 1}, {0.0, 0.5, 0.0}},
      Reading{"x+0.5,y,z", {1, 0, 0, 0, 1, 0, 0, 0, 1}, {0.5, 0.0, 0.0}},
      Reading{"X,Y,Z", {1, 0, 0, 0, 1, 0, 0, 0, 1}, {0.0, 0.0, 0.0}},
      Reading{"x, -y-1/2, z-1/2",
              {1, 0, 0, 0, -1, 0, 0, 0, 1},
              {0.0, -0.5, -0.5}},
  };

  for (auto const &[text, rotation, translation] : readings) {
    INFO(text);
    auto const parsed = from_xyz(text);
    REQUIRE(parsed);
    CHECK(parsed->rotation == rows(rotation));
    CHECK(parsed->translation.isApprox(translation, 1e-12));
  }
}

TEST_CASE("from_xyz rejects text that is not three coordinates", "[xyz]") {
  for (auto const text : std::array<std::string_view, 6>{
           "x,y", "x,y,z,z", "x,y,w", "x,2x,z",
           // Parses, but its rotation has determinant 6: rejected by the
           // unimodularity gate, not by the grammar.
           "-2y+1/2, 3x+1/2, z-y+1/2", ""}) {
    INFO(text);
    CHECK(errored([&] { return from_xyz(text); }));
  }
}

TEST_CASE("to_xyz writes the exact fractions", "[xyz]") {
  auto const written = [](Vector3d const &t) {
    return to_xyz(SymmetryOperation{Matrix3i::Identity(), t});
  };
  CHECK(written({0.0, 0.0, 0.0}) == "x,y,z");
  CHECK(written({1.0 / 3.0, 2.0 / 3.0, 1.0 / 12.0}) == "x+1/3,y+2/3,z+1/12");
  CHECK(written({-0.5, 0.375, 0.25}) == "x-1/2,y+3/8,z+1/4");
  // No exact small denominator: a short decimal instead of a wrong fraction.
  CHECK(written({0.137, 0.0, 0.0}).starts_with("x+0.137"));
}

TEST_CASE("to_xyz writes coefficients and drops the leading plus", "[xyz]") {
  Matrix3i r;
  r << 0, -1, 0, 1, -1, 0, 0, 0, 1;
  CHECK(to_xyz(SymmetryOperation{r, Vector3d(0.0, 0.0, 0.5)}) == "-y,x-y,z+1/2");
  CHECK(to_xyz(SymmetryOperation{Matrix3i::Zero(), Vector3d(0.0, 0.0, 0.0)}) ==
        "0,0,0");
}

TEST_CASE("every database operation round-trips through its triplet",
          "[xyz]") {
  for (int index : std::views::iota(1, kSpaceHallSettings + 1)) {
    auto const hall = test::space_hall(index);
    for (SymmetryOperation const &op : data::operations_from_database(hall)) {
      std::string const text = to_xyz(op);
      INFO(index << ' ' << text);
      auto const back = from_xyz(text);
      REQUIRE(back);
      CHECK(same_operation(*back, op, 1e-9));
    }
  }
}
