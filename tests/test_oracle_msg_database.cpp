// Oracle test for the magnetic space-group database (Phase 7 foundation):
// every UNI number's type metadata and decoded magnetic operations must match
// the reference spglib (spg_get_magnetic_spacegroup_type /
// spg_get_magnetic_symmetry_from_database) bit-for-bit.

#include "oracle.hpp"
#include <seitz/core/operation_set.hpp>

#include <seitz/data/msg_database.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <string>
#include <vector>

namespace {

using seitz::MagneticOperations;
using seitz::Matrix3i;
using seitz::Vector3d;
using seitz::test::space_hall;
using seitz::test::uni_number;

// Reference magnetic operations for (uni_number, hall_number) via spglib's
// spg_get_magnetic_symmetry_from_database (max 384 operations).
MagneticOperations reference_magnetic_operations(int uni_number,
                                                 int hall_number) {
  std::array<int, 384 * 9> rot{};
  std::array<double, 384 * 3> trans{};
  std::array<int, 384> timerev{};
  int const n = spg_get_magnetic_symmetry_from_database(
      reinterpret_cast<int(*)[3][3]>(rot.data()),
      reinterpret_cast<double(*)[3]>(trans.data()), timerev.data(), uni_number,
      hall_number);
  std::vector<seitz::MagneticSymmetryOperation> ops;
  ops.reserve(static_cast<std::size_t>(n));
  for (int s = 0; s < n; ++s) {
    Matrix3i r;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        r(i, j) = rot[static_cast<std::size_t>(9 * s + 3 * i + j)];
    Vector3d t(trans[static_cast<std::size_t>(3 * s + 0)],
               trans[static_cast<std::size_t>(3 * s + 1)],
               trans[static_cast<std::size_t>(3 * s + 2)]);
    ops.push_back({{r, t}, timerev[static_cast<std::size_t>(s)] != 0});
  }
  return MagneticOperations{std::move(ops)};
}

void require_same_operations(MagneticOperations const &got,
                             MagneticOperations const &ref) {
  REQUIRE(got.size() == ref.size());
  for (std::size_t i = 0; i < ref.size(); ++i) {
    INFO("operation index " << i);
    REQUIRE(got[i].spatial.rotation == ref[i].spatial.rotation);
    REQUIRE((got[i].spatial.translation - ref[i].spatial.translation)
                .cwiseAbs()
                .maxCoeff() < 1e-12);
    REQUIRE(got[i].time_reversal == ref[i].time_reversal);
  }
}

} // namespace

TEST_CASE("magnetic_spacegroup_type matches reference for all UNI numbers",
          "[oracle][msg_database]") {
  for (int uni = 1; uni <= seitz::data::kNumUniNumbers; ++uni) {
    INFO("UNI number " << uni);
    auto const &got =
        seitz::data::magnetic_spacegroup_type(uni_number(uni));
    SpglibMagneticSpacegroupType const ref =
        spg_get_magnetic_spacegroup_type(uni);

    REQUIRE(got.uni_number == ref.uni_number);
    REQUIRE(got.litvin_number == ref.litvin_number);
    REQUIRE(got.bns_number == std::string(ref.bns_number));
    REQUIRE(got.og_number == std::string(ref.og_number));
    REQUIRE(got.number == ref.number);
    REQUIRE(got.type == ref.type);
  }
}

TEST_CASE("magnetic operations match reference (default Hall setting)",
          "[oracle][msg_database]") {
  for (int uni = 1; uni <= seitz::data::kNumUniNumbers; ++uni) {
    INFO("UNI number " << uni);
    auto const &got =
        seitz::data::magnetic_operations_from_database(uni_number(uni));
    auto const ref = reference_magnetic_operations(uni, 0);
    REQUIRE(!got.empty());
    require_same_operations(got, ref);
  }
}

TEST_CASE("magnetic operations match reference across all Hall settings",
          "[oracle][msg_database]") {
  for (int index = 1; index <= seitz::kSpaceHallSettings; ++index) {
    auto const [first, last] =
        seitz::data::uni_candidates(space_hall(index));
    for (int uni = first.value(); uni <= last.value(); ++uni) {
      INFO("UNI number " << uni << ", Hall number " << index);
      auto const &got = seitz::data::magnetic_operations_from_database(
          uni_number(uni), space_hall(index));
      auto const ref = reference_magnetic_operations(uni, index);
      require_same_operations(got, ref);
    }
  }
}
