#include <seitz/core/cell.hpp>

#include "helpers.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <vector>

using namespace seitz;

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
    CHECK(position.isApprox(
        c.position(static_cast<Index>(seen_types.size()) - 1)));
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

TEST_CASE("Cell::supercell repeats the atoms over the new lattice points",
          "[cell]") {
  Cell const motif = test::rocksalt_motif(4.0);
  Cell const sc = test::must(motif.supercell(Vector3i(2, 2, 2)));
  CHECK(sc.size() == 16);
  CHECK(sc.lattice().volume() == Catch::Approx(8.0 * motif.lattice().volume()));
  CHECK(std::ranges::count(sc.types(), 0) == 8);
  CHECK((sc.positions().array() >= 0.0).all());
  CHECK((sc.positions().array() < 1.0).all());

  // The same space group as the motif, and the motif back as the primitive.
  auto const motif_hall = test::must(test::dataset_of(motif)).hall;
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(sc);
  CHECK(test::must(analyzer.hall()) == motif_hall);
  CHECK(test::must(analyzer.primitive_cell()).size() == motif.size());
}

TEST_CASE("Cell::transformed takes a non-diagonal unimodular basis", "[cell]") {
  Cell const motif = test::rocksalt_motif(4.0);
  Matrix3i basis;
  basis << 1, 1, 0, -1, 1, 0, 0, 0, 1; // a' = a - b, b' = a + b
  Cell const cell = test::must(motif.transformed(basis, Vector3d(0.5, 0, 0)));
  CHECK(cell.size() == 4);
  CHECK(cell.lattice().matrix().isApprox(motif.lattice().matrix() *
                                         basis.cast<double>()));
  CHECK(test::must(test::dataset_of(cell)).hall ==
        test::must(test::dataset_of(motif)).hall);
}

TEST_CASE("Cell::translated keeps the symmetry and folds into the cell",
          "[cell]") {
  Cell const motif = test::rocksalt_motif(4.0);
  Cell const moved = motif.translated(Vector3d(0.7, -0.2, 1.3));
  CHECK(moved.size() == motif.size());
  CHECK((moved.positions().array() >= 0.0).all());
  CHECK((moved.positions().array() < 1.0).all());
  CHECK(moved.position(0).isApprox(Vector3d(0.7, 0.8, 0.3)));
  CHECK(test::must(test::dataset_of(moved)).hall ==
        test::must(test::dataset_of(motif)).hall);
}

TEST_CASE("Cell::transformed rejects singular and axis-mixing bases",
          "[cell]") {
  Cell const motif = test::rocksalt_motif(4.0);
  CHECK(test::errored([&] { return motif.transformed(Matrix3i::Zero()); }));

  // A layer: the aperiodic c may be kept or flipped, never mixed in.
  Cell const layer = motif.with_periodicity(aperiodic_along(2));
  Matrix3i mixing;
  mixing << 1, 0, 1, 0, 1, 0, 0, 0, 1;
  CHECK(test::errored([&] { return layer.transformed(mixing); }));
  Matrix3i flipped = Matrix3i::Identity();
  flipped(2, 2) = -1;
  Cell const kept = test::must(layer.transformed(flipped));
  CHECK(kept.periodicity() == layer.periodicity());
  // The aperiodic coordinate is negated, not folded.
  CHECK(kept.position(1)[2] == Catch::Approx(-0.5));
  Cell const doubled =
      test::must(layer.transformed(Vector3i(2, 1, 1).asDiagonal()));
  CHECK(doubled.size() == 4);
}
