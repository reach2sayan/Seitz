#pragma once

#include <seitz/analysis/magnetic_symmetry_analyzer.hpp>
#include <seitz/analysis/symmetry_analyzer.hpp>
#include <seitz/core/error.hpp>
#include <seitz/core/fractional.hpp>
#include <seitz/core/keys.hpp>
#include <seitz/core/operation_set.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/group/space_group.hpp>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <ranges>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

// Shared Result<T> plumbing for the test suites, so each file does not carry
// its own copy of the same two try_handle_all wrappers.
namespace seitz::test {

// Shorthand for the validated keys: the tests name literal settings by the
// dozen, and `*HallNumber::of(GroupFamily::space, 1)` at each of them buries
// what is being tested.
[[nodiscard]] constexpr HallNumber space_hall(int index) noexcept {
  return *HallNumber::of(GroupFamily::space, index);
}
[[nodiscard]] constexpr HallNumber layer_hall(int index) noexcept {
  return *HallNumber::of(GroupFamily::layer, index);
}
[[nodiscard]] constexpr UniNumber uni_number(int number) noexcept {
  return *UniNumber::of(number);
}

// The determination of one cell, the way callers reach it now. A shorthand
// only: SymmetryAnalyzer::from_cell(...).dataset() at two dozen call sites
// would bury what each test is actually checking.
// The analyzer's accessors hand back references into its memo, so the analyzer
// has to outlive them: these own one for the duration and copy out.
[[nodiscard]] inline Result<analysis::Dataset>
dataset_of(Cell const &cell, Tolerance const &tol = {}) {
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(cell, tol);
  BOOST_LEAF_AUTO(ds, analyzer.dataset());
  return ds;
}

[[nodiscard]] inline Result<analysis::MagneticDataset>
magnetic_dataset_of(MagneticCell const &cell,
                    MagneticTolerance const &tol = {}) {
  auto const analyzer =
      analysis::MagneticSymmetryAnalyzer::from_cell(cell, tol);
  BOOST_LEAF_AUTO(ds, analyzer.dataset());
  return ds;
}

// The success value of `r`, failing the current test if it carries an error.
// Always by value: the analyzers hand back references into their memo, and
// returning T unchanged would bind a reference to try_handle_all's temporary.
template <class T> std::remove_cvref_t<T> must(Result<T> r) {
  using Value = std::remove_cvref_t<T>;
  return leaf::try_handle_all([&]() -> Result<Value> { return std::move(r); },
                              [](leaf::error_info const &) -> Value {
                                FAIL("unexpected error result");
                                throw std::logic_error("unreachable");
                              });
}

// Whether `make()` produced an error — the expected-failure counterpart of
// must(). Takes the call rather than the Result so a BOOST_LEAF_AUTO chain can
// be written inline.
template <ResultProducer F> [[nodiscard]] bool errored(F &&make) {
  return leaf::try_handle_all(
      [&]() -> Result<bool> {
        if (auto r = std::forward<F>(make)(); !r) {
          return r.error();
        }
        return false;
      },
      [](leaf::error_info const &) { return true; });
}

// The space (layer) group of international number `n` in its default setting,
// failing the test if the lookup errors. Flyweights: the pointer is stable.
[[nodiscard]] inline group::SpaceGroup const *space_group(int n) {
  return must(group::SpaceGroup::from_number(GroupFamily::space, n));
}
[[nodiscard]] inline group::SpaceGroup const *layer_group(int n) {
  return must(group::SpaceGroup::from_number(GroupFamily::layer, n));
}

// A cell from literal fractional rows. The cubic-edge overload is the common
// shape of the hand-built magnetic and spin structures.
[[nodiscard]] inline Cell make_cell(Matrix3d const &lattice,
                                    std::vector<std::array<double, 3>> const &pos,
                                    Types types) {
  Positions p(static_cast<Index>(pos.size()), 3);
  for (auto const [i, row] : pos | std::views::enumerate) {
    p.row(i) = Eigen::Map<Eigen::RowVector3d const>(row.data());
  }
  return Cell(Lattice{lattice}, p, std::move(types));
}
[[nodiscard]] inline Cell make_cell(double a,
                                    std::vector<std::array<double, 3>> const &pos,
                                    Types types) {
  return make_cell(Matrix3d{a * Matrix3d::Identity()}, pos, std::move(types));
}

// The two-atom rocksalt motif (CsCl arrangement) in a cubic cell of edge `a`:
// the seed the supercell-shaped tests and benchmarks grow from.
[[nodiscard]] inline Cell rocksalt_motif(double a, Types types = {0, 1}) {
  return make_cell(a, {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}}, std::move(types));
}

// Whether two fractional points coincide modulo 1 within `tol` (inclusive).
[[nodiscard]] inline bool coincide_mod1(Vector3d const &a, Vector3d const &b,
                                        double tol) {
  return math::nearest_offset(Vector3d(a - b)).cwiseAbs().maxCoeff() <= tol;
}

// Equal as sets: rotations exact, translations modulo the lattice. Reference
// implementations order their operations differently than we discover them.
[[nodiscard]] inline bool same_operation_set(auto const &a, auto const &b,
                                             double tol) {
  return std::ranges::size(a) == std::ranges::size(b) &&
         std::ranges::all_of(a, [&](SymmetryOperation const &x) {
           return std::ranges::any_of(b, [&](SymmetryOperation const &y) {
             return same_operation(x, y, tol);
           });
         });
}

// Two lattices agree up to a rigid rotation iff their metric tensors do.
[[nodiscard]] inline bool same_metric(Matrix3d const &a, Matrix3d const &b,
                                      double tol) {
  return (Lattice{a}.metric() - Lattice{b}.metric()).cwiseAbs().maxCoeff() <
         tol;
}

// Two bases span the same lattice iff a^-1 b is an integer unimodular matrix.
[[nodiscard]] inline bool same_lattice(Matrix3d const &a, Matrix3d const &b,
                                       double tol = 1e-6) {
  Matrix3d const rel = a.inverse() * b;
  Matrix3i const rounded = math::round_to_int(rel);
  return (rel - rounded.cast<double>()).cwiseAbs().maxCoeff() <= tol &&
         std::abs(rounded.determinant()) == 1;
}

// Every atom displaced by a random Cartesian vector of amplitude
// `cart_amplitude` per component.
[[nodiscard]] inline Cell jitter(Cell const &cell, std::mt19937 &rng,
                                 double cart_amplitude) {
  Matrix3d const lattice_inv = cell.lattice().matrix().inverse();
  std::uniform_real_distribution<double> dist(-cart_amplitude, cart_amplitude);
  Positions positions = cell.positions();
  for (auto row : positions.rowwise()) {
    row += (lattice_inv * Vector3d::NullaryExpr([&](Index) { return dist(rng); }))
               .transpose();
  }
  return Cell(cell.lattice(), positions, cell.types(), cell.periodicity());
}

// The k x k x k rocksalt supercell (edge 4k angstrom, 2 k^3 atoms), jittered by
// `noise` in fractional units of the supercell.
[[nodiscard]] inline Cell rocksalt_supercell(int k, double noise,
                                             std::mt19937 &rng) {
  return jitter(must(rocksalt_motif(4.0).supercell(Vector3i::Constant(k))), rng,
                noise * 4.0 * k);
}

} // namespace seitz::test
