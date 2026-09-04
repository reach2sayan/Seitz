// Phase 9 hardening: degenerate inputs must produce a clean Result error (the
// right tag) rather than UB. No reference spglib needed — this exercises the
// port's own validation/guards.

#include "symmetry/search.hpp"
#include <cppcrystal/analysis/magnetic_symmetry_analyzer.hpp>
#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>

#include "helpers.hpp"

#include <catch2/catch_test_macros.hpp>

#include <vector>

namespace {

using namespace cppcrystal;

enum class Err { none, empty, invalid_lattice, other };

// Run `f` (returning a Result) inside a LEAF handling scope and report which
// error tag, if any, it produced.
template <class F> Err classify(F f) {
  Err err = Err::none;
  leaf::try_handle_all(
      [&]() -> leaf::result<void> {
        auto r = f();
        if (!r) {
          return r.error();
        }
        return {};
      },
      [&](e_empty_cell const &) { err = Err::empty; },
      [&](e_invalid_lattice const &) { err = Err::invalid_lattice; },
      [&]() { err = Err::other; });
  return err;
}

Cell empty_cell() {
  return Cell(Lattice{4.0 * Matrix3d::Identity()}, Positions(0, 3), Types{});
}

// Two collinear basis vectors -> determinant 0 (a degenerate "lattice").
Cell singular_cell() {
  Matrix3d lattice;
  lattice.col(0) = Vector3d(4, 0, 0);
  lattice.col(1) = Vector3d(4, 0, 0); // == col(0)
  lattice.col(2) = Vector3d(0, 0, 4);
  Positions p(1, 3);
  p.row(0) = Eigen::RowVector3d(0, 0, 0);
  return Cell(Lattice{lattice}, p, Types{0});
}

Cell valid_cell() {
  Positions p(1, 3);
  p.row(0) = Eigen::RowVector3d(0, 0, 0);
  return Cell(Lattice{4.0 * Matrix3d::Identity()}, p, Types{0});
}

} // namespace

TEST_CASE("empty cell is rejected with e_empty_cell", "[error]") {
  CHECK(classify([] { return test::dataset_of(empty_cell(), {1e-5}); }) ==
        Err::empty);
  CHECK(classify([] {
          Cell const cell = empty_cell();
          return symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}
              .operations();
        }) == Err::empty);
  MagneticCell const m(empty_cell(), SiteTensors{CollinearTensors{}},
                       TensorKind::axial);
  CHECK(classify([&] { return test::magnetic_dataset_of(m, {1e-5}); }) ==
        Err::empty);
}

TEST_CASE("singular lattice is rejected with e_invalid_lattice", "[error]") {
  CHECK(classify([] { return test::dataset_of(singular_cell(), {1e-5}); }) ==
        Err::invalid_lattice);
  CHECK(classify([] {
          Cell const cell = singular_cell();
          return symmetry::SymmetrySearch<GroupFamily::space>{cell, {1e-5}}
              .operations();
        }) == Err::invalid_lattice);
  MagneticCell const m(singular_cell(), SiteTensors{CollinearTensors{1.0}},
                       TensorKind::axial);
  CHECK(classify([&] { return test::magnetic_dataset_of(m, {1e-5}); }) ==
        Err::invalid_lattice);
}

// A bad mesh is now unrepresentable rather than something every grid function
// has to guard: Mesh::of rejects it once, at the only entry point. Covered by
// the static_asserts in kpoint/mesh.hpp and by test_kpoint.cpp.
