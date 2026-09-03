#include <cppcrystal/core/cell.hpp>

#include <catch2/catch_test_macros.hpp>

#include <vector>

using namespace cppcrystal;

TEST_CASE("Cell exposes size, types, positions", "[cell]") {
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  Cell c(Lattice{Matrix3d::Identity() * 3.0}, pos, {0, 1});
  CHECK(c.size() == 2);
  CHECK(c.type(1) == 1);
  CHECK(c.position(1).isApprox(Vector3d(0.5, 0.5, 0.5)));
  CHECK(c.periodicity() == all_periodic());
  CHECK(c.lattice().matrix().isApprox(Matrix3d::Identity() * 3.0));
}

TEST_CASE("Cell atoms() pairs each position with its type", "[cell]") {
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  Cell const c(Lattice{Matrix3d::Identity() * 3.0}, pos, {7, 9});

  std::vector<int> seen_types;
  for (auto const &[position, type] : c.atoms()) {
    seen_types.push_back(type);
    CHECK(position.isApprox(c.position(static_cast<Index>(seen_types.size()) - 1)));
  }
  CHECK(seen_types == std::vector<int>{7, 9});
}

TEST_CASE("Cell with_* builders replace one facet and keep the atoms",
          "[cell]") {
  Positions pos(1, 3);
  pos.row(0) << 0.25, 0.25, 0.25;
  Cell const c(Lattice{Matrix3d::Identity() * 3.0}, pos, {0});

  Cell const relatticed = c.with_lattice(Lattice{Matrix3d::Identity() * 5.0});
  CHECK(relatticed.lattice().volume() == 125.0);
  CHECK(relatticed.positions().isApprox(c.positions()));
  CHECK(relatticed.types() == c.types());

  CellPeriodicity const layer{AxisKind::periodic, AxisKind::periodic,
                              AxisKind::aperiodic};
  Cell const layered = c.with_periodicity(layer);
  CHECK(layered.periodicity() == layer);
  CHECK(aperiodic_axis(layered.periodicity()) == 2);
  CHECK(layered.lattice().matrix().isApprox(c.lattice().matrix()));
}
