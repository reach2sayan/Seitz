#include <spglib/data/spg_database.hpp>

#include <catch2/catch_test_macros.hpp>

using namespace spglib;

TEST_CASE("decode_operation recovers identity and inversion", "[database]") {
  // Hall 1 is P1: a single identity operation.
  auto p1 = data::operations_from_database(1);
  REQUIRE(p1.size() == 1);
  CHECK(p1[0].rotation == Matrix3i::Identity());
  CHECK(p1[0].translation.cwiseAbs().maxCoeff() < 1e-12);

  // Hall 2 is P-1: identity + inversion.
  auto p1bar = data::operations_from_database(2);
  REQUIRE(p1bar.size() == 2);
  bool has_inversion = false;
  for (auto const &op : p1bar)
    if (op.rotation == -Matrix3i::Identity())
      has_inversion = true;
  CHECK(has_inversion);
}

TEST_CASE("spacegroup_type returns trimmed symbols", "[database]") {
  auto t1 = data::spacegroup_type(1);
  CHECK(t1.number == 1);
  CHECK(t1.international_short == "P1");
  CHECK(t1.centering == data::Centering::primitive);

  auto t_last = data::spacegroup_type(530);
  CHECK(t_last.number == 230); // hall 530 is the last setting of Ia-3d
}

TEST_CASE("out-of-range hall numbers are handled", "[database]") {
  CHECK(data::operations_from_database(0).empty());
  CHECK(data::operations_from_database(531).empty());
  CHECK(data::spacegroup_type(0).number == 0);
}
