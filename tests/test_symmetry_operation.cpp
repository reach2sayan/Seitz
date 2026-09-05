#include "core/matrix_order.hpp"
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <iterator>

using namespace seitz;
using Catch::Approx;

namespace {
SymmetryOperation rot_z_90() {
  Matrix3i r;
  r << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  return {r, Vector3d(0.0, 0.0, 0.0)};
}
} // namespace

TEST_CASE("apply computes rot . x + trans", "[symop]") {
  SymmetryOperation op{Matrix3i::Identity(), Vector3d(0.5, 0.0, 0.0)};
  CHECK(op.apply(Vector3d(0.1, 0.2, 0.3)).isApprox(Vector3d(0.6, 0.2, 0.3)));
}

TEST_CASE("composition matches sequential application", "[symop]") {
  SymmetryOperation a = rot_z_90();
  a.translation = Vector3d(0.5, 0.0, 0.0);
  SymmetryOperation b{Matrix3i::Identity(), Vector3d(0.0, 0.25, 0.0)};
  SymmetryOperation ab = a * b;
  Vector3d x(0.1, 0.2, 0.3);
  CHECK(ab.apply(x).isApprox(a.apply(b.apply(x))));
}

TEST_CASE("inverse undoes the operation", "[symop]") {
  SymmetryOperation a = rot_z_90();
  a.translation = Vector3d(0.5, 0.1, 0.0);
  auto inv = a.inverse();
  REQUIRE(inv.has_value());
  SymmetryOperation id = a * *inv;
  CHECK(id.rotation == Matrix3i::Identity());
  CHECK(id.translation.cwiseAbs().maxCoeff() == Approx(0.0).margin(1e-12));
}

TEST_CASE("same_operation respects symprec on translations", "[symop]") {
  SymmetryOperation a = rot_z_90();
  a.translation = Vector3d(0.5, 0.0, 0.0);
  SymmetryOperation b = a;
  b.translation += Vector3d(1.0, -2.0, 3.0); // differ by a lattice vector
  CHECK(same_operation(a, b, 1e-5));
  b.translation = a.translation + Vector3d(0.01, 0.0, 0.0);
  CHECK_FALSE(same_operation(a, b, 1e-5));
}

TEST_CASE("index_by_rotation keeps original order within a rotation",
          "[symop][matrix_order]") {
  SymmetryOperation const id{Matrix3i::Identity(), Vector3d::Zero()};
  SymmetryOperation const id_shifted{Matrix3i::Identity(),
                                     Vector3d(0.5, 0.5, 0.5)};
  Operations const ops{
      std::vector<SymmetryOperation>{rot_z_90(), id, id_shifted}};

  auto const by_rot = index_by_rotation(ops, &SymmetryOperation::rotation);
  REQUIRE(by_rot.size() == 3);
  CHECK(by_rot.find(rot_z_90().rotation)->second == 0);
  auto const [lo, hi] = by_rot.equal_range(Matrix3i::Identity());
  REQUIRE(std::distance(lo, hi) == 2);
  CHECK(lo->second == 1);
  CHECK(std::next(lo)->second == 2);
  CHECK(by_rot.find(-Matrix3i::Identity()) == by_rot.end());

  CHECK(has_duplicate_rotation(ops, &SymmetryOperation::rotation));
  CHECK_FALSE(has_duplicate_rotation(
      Operations{std::vector<SymmetryOperation>{rot_z_90(), id}},
      &SymmetryOperation::rotation));
  CHECK(rotation_set(ops, &SymmetryOperation::rotation).size() == 2);
}
