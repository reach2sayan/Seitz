#include <cppcrystal/data/sitesym_database.hpp>
#include <cppcrystal/data/sitesym_tables.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace cppcrystal;
using namespace cppcrystal::data;
using cppcrystal::test::space_hall;

TEST_CASE("wyckoff_coordinate decodes the documented entries", "[sitesym]") {
  // Index 1 (P1 general position): identity rotation, zero translation.
  auto const w1 = wyckoff_coordinate(1);
  CHECK(w1.rotation == Matrix3i::Identity());
  CHECK(w1.translation.isZero());
  CHECK(w1.multiplicity == 1);

  // Index 3 (a P-1 inversion centre): the rotation is the zero matrix (a fully
  // fixed point) at (1/2, 1/2, 1/2).
  auto const w3 = wyckoff_coordinate(3);
  CHECK(w3.rotation == Matrix3i::Zero());
  CHECK(w3.translation.isApprox(Vector3d(0.5, 0.5, 0.5)));
  CHECK(w3.multiplicity == 1);
}

TEST_CASE("wyckoff_coordinate rotation/translation stay in encodable ranges",
          "[sitesym]") {
  // First column entries in {-2..2}, other entries in {-1..1}; translations are
  // multiples of 1/24 in [0, 1). Spot-check across the whole 3D table.
  int const last = kPositionWyckoff[530]; // first index of Hall 530's range end
  for (int i = 1; i < last; ++i) {
    auto const w = wyckoff_coordinate(i);
    for (int r = 0; r < 3; ++r) {
      CHECK(w.rotation(r, 0) >= -2);
      CHECK(w.rotation(r, 0) <= 2);
      CHECK(w.rotation(r, 1) >= -1);
      CHECK(w.rotation(r, 1) <= 1);
      CHECK(w.rotation(r, 2) >= -1);
      CHECK(w.rotation(r, 2) <= 1);
    }
    for (int k = 0; k < 3; ++k) {
      double const x = w.translation[k] * 24.0;
      CHECK(x >= -0.001);
      CHECK(x < 24.001);
      CHECK(std::abs(x - std::round(x)) < 1e-9); // a multiple of 1/24
    }
    CHECK(w.multiplicity > 0);
  }
}

TEST_CASE("wyckoff_indices match the documented Wyckoff counts", "[sitesym]") {
  // P1 (Hall 1): 1 Wyckoff position starting at global index 1.
  auto const r1 = wyckoff_indices(space_hall(1));
  CHECK(r1.start == 1);
  CHECK(r1.count == 1);

  // P-1 (Hall 2): 9 Wyckoff positions.
  auto const r2 = wyckoff_indices(space_hall(2));
  CHECK(r2.start == 2);
  CHECK(r2.count == 9);

  // Every Hall number's range is non-negative, ordered and within the table.
  int const table_size = static_cast<int>(kCoordinatesFirst.size());
  for (int index = 1; index <= 530; ++index) {
    auto const r = wyckoff_indices(space_hall(index));
    INFO("hall " << index);
    CHECK(r.count > 0);
    CHECK(r.start >= 1);
    CHECK(r.start + r.count <= table_size);
  }
}

TEST_CASE("site_symmetry_symbol returns the pre-trimmed symbol", "[sitesym]") {
  CHECK(site_symmetry_symbol(1) == "1");
  CHECK(site_symmetry_symbol(3) == "-1");
  CHECK(site_symmetry_symbol(0).empty()); // empty dummy entry
}
