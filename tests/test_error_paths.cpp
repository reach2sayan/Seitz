// Phase 9 hardening: degenerate inputs must produce a clean Result error (the
// right tag) rather than UB. No reference spglib needed — this exercises the
// port's own validation/guards.

#include <cppcrystal/core/cell.hpp>
#include <cppcrystal/core/error.hpp>
#include <cppcrystal/core/magnetic_cell.hpp>
#include <cppcrystal/dataset.hpp>
#include <cppcrystal/kpoint/grid.hpp>
#include <cppcrystal/kpoint/reciprocal_mesh.hpp>
#include <cppcrystal/magnetic_dataset.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>

#include <catch2/catch_test_macros.hpp>

#include <vector>

namespace {

using namespace cppcrystal;

enum class Err { none, empty, invalid_lattice, invalid_mesh, other };

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
      [&](e_invalid_mesh const &) { err = Err::invalid_mesh; },
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
  CHECK(classify([] { return get_dataset(empty_cell(), {1e-5}); }) == Err::empty);
  CHECK(classify([] { return symmetry::find_symmetry(empty_cell(), {1e-5}); }) ==
        Err::empty);
  MagneticCell const m(empty_cell(), SiteTensors{CollinearTensors{}},
                       TensorKind::axial);
  CHECK(classify([&] { return get_magnetic_dataset(m, {1e-5}); }) ==
        Err::empty);
}

TEST_CASE("singular lattice is rejected with e_invalid_lattice", "[error]") {
  CHECK(classify([] { return get_dataset(singular_cell(), {1e-5}); }) ==
        Err::invalid_lattice);
  CHECK(classify([] {
          return symmetry::find_symmetry(singular_cell(), {1e-5});
        }) == Err::invalid_lattice);
  MagneticCell const m(singular_cell(), SiteTensors{CollinearTensors{1.0}},
                       TensorKind::axial);
  CHECK(classify([&] { return get_magnetic_dataset(m, {1e-5}); }) ==
        Err::invalid_lattice);
}

TEST_CASE("non-positive mesh is rejected with e_invalid_mesh", "[error]") {
  CHECK(classify([] {
          return kpoint::ir_reciprocal_mesh(valid_cell(), Vector3i(4, 0, 4),
                                            Vector3i::Zero(), TimeReversal::on);
        }) == Err::invalid_mesh);
  CHECK(classify([] {
          return kpoint::ir_reciprocal_mesh(valid_cell(), Vector3i(4, -2, 4),
                                            Vector3i::Zero(), TimeReversal::on);
        }) == Err::invalid_mesh);
}

TEST_CASE("k-point value functions degrade instead of UB on a bad mesh",
          "[error]") {
  // Zero / negative mesh must not divide-by-zero or over-allocate.
  CHECK(kpoint::all_grid_addresses(Vector3i(4, 0, 4)).empty());
  CHECK(kpoint::all_grid_addresses(Vector3i(-1, 4, 4)).empty());
  CHECK(kpoint::grid_point_from_address(Vector3i(1, 1, 1),
                                        Vector3i(0, 0, 0)) == 0U);
  // A valid mesh still works (sanity).
  CHECK(kpoint::all_grid_addresses(Vector3i(2, 2, 2)).size() == 8U);
}
