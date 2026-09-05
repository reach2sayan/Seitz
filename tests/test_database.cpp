#include <seitz/data/element_data.hpp>
#include <seitz/data/msg_database.hpp>
#include <seitz/data/spg_database.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

using namespace seitz;
using seitz::test::space_hall;
using seitz::test::uni_number;

TEST_CASE("operations_from_database returns identity and inversion",
          "[database]") {
  // Hall 1 is P1: a single identity operation.
  auto const &p1 = data::operations_from_database(space_hall(1));
  REQUIRE(p1.size() == 1);
  CHECK(p1[0].rotation == Matrix3i::Identity());
  CHECK(p1[0].translation.cwiseAbs().maxCoeff() < 1e-12);

  // Hall 2 is P-1: identity + inversion.
  auto const &p1bar = data::operations_from_database(space_hall(2));
  REQUIRE(p1bar.size() == 2);
  bool has_inversion = false;
  for (auto const &op : p1bar)
    if (op.rotation == -Matrix3i::Identity())
      has_inversion = true;
  CHECK(has_inversion);
}

TEST_CASE("spacegroup_type returns trimmed symbols", "[database]") {
  auto const &t1 = data::spacegroup_type(space_hall(1));
  CHECK(t1.number == 1);
  CHECK(t1.international_short == "P1");
  CHECK(t1.centering == data::Centering::primitive);

  auto const &t_last = data::spacegroup_type(space_hall(530));
  CHECK(t_last.number == 230); // hall 530 is the last setting of Ia-3d
}

TEST_CASE("a Hall number cannot name a setting that does not exist",
          "[database]") {
  // The old sentinel-row fallback is gone: out-of-range is now unrepresentable
  // rather than something the catalog has to answer for.
  CHECK_FALSE(HallNumber::of(GroupFamily::space, 0).has_value());
  CHECK_FALSE(HallNumber::of(GroupFamily::space, 531).has_value());
  CHECK_FALSE(HallNumber::of(GroupFamily::space, -1).has_value());
  CHECK_FALSE(HallNumber::of(GroupFamily::layer, 117).has_value());
  CHECK(HallNumber::of(GroupFamily::layer, 116).has_value());
  CHECK_FALSE(UniNumber::of(0).has_value());
  CHECK_FALSE(UniNumber::of(1652).has_value());
}

TEST_CASE("atomic_number inverts element_symbol for every element",
          "[database]") {
  for (int z = 1; z <= data::kNumElements; ++z) {
    CHECK(data::atomic_number(*data::element_symbol(z)) == z);
  }
  CHECK_FALSE(data::atomic_number("Xx"));
  CHECK_FALSE(data::atomic_number(""));
}

TEST_CASE("a (uni, hall) pairing outside the uni's settings yields nothing",
          "[database]") {
  // UNI 1 (P1) has a single Hall setting (1); any other Hall number is not
  // one of its settings, and the query must not read past that row.
  CHECK_FALSE(
      data::magnetic_operations_from_database(uni_number(1), space_hall(1))
          .empty());
  CHECK(data::magnetic_operations_from_database(uni_number(1), space_hall(2))
            .empty());
  CHECK(data::magnetic_operations_from_database(uni_number(1), space_hall(530))
            .empty());
  CHECK(
      data::magnetic_std_transformations(uni_number(1), space_hall(2)).empty());
}
