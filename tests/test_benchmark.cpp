// End-to-end timing for the three user-visible workloads: symmetry
// determination, reciprocal-mesh reduction, and structure generation. All of it
// is tagged [!benchmark], so it is skipped by a plain `ctest` run and only
// executes under an explicit tag filter:
//
//   ./cppcrystal_tests "[!benchmark]"
//
// The cells are built here rather than loaded from the reference corpus: the
// corpus loader needs SPGLIB_REF_DATA_DIR, which only the oracle target
// defines, and a benchmark that cannot run in the default build is a benchmark
// nobody runs.

#include <cppcrystal/analysis/symmetry_analyzer.hpp>
#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/generate/generator.hpp>
#include <cppcrystal/group/space_group.hpp>
#include <cppcrystal/kpoint/mesh.hpp>

#include "helpers.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <random>
#include <vector>

using namespace cppcrystal;

namespace {

// A k x k x k supercell of a 2-atom rocksalt-like motif, with the second
// species split across two types so that type filtering matters. Same shape as
// the helper in test_overlap.cpp; duplicated rather than shared because the two
// files exercise it for different reasons and a shared header would tie the
// benchmark's cell to the overlap test's.
Cell supercell(int k, double noise, std::mt19937 &rng) {
  Index const n = 2 * k * k * k;
  Matrix3d const lattice = Matrix3d::Identity() * (4.0 * k);
  Positions pos(n, 3);
  Types types;
  types.reserve(static_cast<std::size_t>(n));
  std::uniform_real_distribution<double> jitter(-noise, noise);
  Index row = 0;
  for (int a = 0; a < k; ++a) {
    for (int b = 0; b < k; ++b) {
      for (int c = 0; c < k; ++c) {
        Vector3d const base(a, b, c);
        pos.row(row++) = (base / k).transpose();
        types.push_back(0);
        pos.row(row++) = ((base + Vector3d(0.5, 0.5, 0.5)) / k).transpose();
        types.push_back((a + b + c) % 2);
      }
    }
  }
  for (Index i = 0; i < n; ++i) {
    for (Index d = 0; d < 3; ++d) {
      pos(i, d) += jitter(rng);
    }
  }
  return Cell(Lattice{lattice}, pos, types);
}

// A low-symmetry cell of the same size, so the timings cover both the fast path
// (a highly symmetric cell, where the Hall search hits early) and the slow one
// (few operations survive, so every candidate setting is tried).
Cell triclinic_cell(int n_atoms, std::mt19937 &rng) {
  Matrix3d lattice;
  lattice.col(0) << 5.1, 0.0, 0.0;
  lattice.col(1) << 1.3, 6.4, 0.0;
  lattice.col(2) << 0.7, 1.1, 7.9;
  Positions pos(n_atoms, 3);
  Types types;
  types.reserve(static_cast<std::size_t>(n_atoms));
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  for (Index i = 0; i < n_atoms; ++i) {
    pos.row(i) << unit(rng), unit(rng), unit(rng);
    types.push_back(static_cast<int>(i) % 3);
  }
  return Cell(Lattice{lattice}, pos, types);
}

// The rotations a determination yields, as the reciprocal mesh wants them.
std::vector<Matrix3i> rotations_of(Cell const &cell) {
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(cell);
  return test::must(analyzer.operations()).rotations();
}

// One determination, timed end to end. The analyzer is built inside so each
// iteration pays for the whole pipeline rather than reading a warm memo; it is
// named rather than a temporary because dataset() is &-qualified.
bool determine(Cell const &cell) {
  auto const analyzer = analysis::SymmetryAnalyzer::from_cell(cell);
  return analyzer.dataset().has_value();
}

} // namespace

TEST_CASE("determine() over supercells", "[!benchmark]") {
  std::mt19937 rng(7);
  Cell const small = supercell(2, 1e-4, rng);   // 16 atoms
  Cell const medium = supercell(4, 1e-4, rng);  // 128 atoms
  Cell const large = supercell(6, 1e-4, rng);   // 432 atoms

  BENCHMARK("cubic supercell, 16 atoms") {
    return determine(small);
  };
  BENCHMARK("cubic supercell, 128 atoms") {
    return determine(medium);
  };
  BENCHMARK("cubic supercell, 432 atoms") {
    return determine(large);
  };
}

TEST_CASE("determine() on a low-symmetry cell", "[!benchmark]") {
  std::mt19937 rng(11);
  // The P1 path: the lattice point group is trivial, so the Hall search walks
  // the candidate list instead of hitting early.
  Cell const triclinic = triclinic_cell(48, rng);

  BENCHMARK("triclinic, 48 atoms") {
    return determine(triclinic);
  };
}

TEST_CASE("reciprocal-mesh reduction", "[!benchmark]") {
  std::mt19937 rng(13);
  Cell const cell = supercell(2, 0.0, rng);
  std::vector<Matrix3i> const rotations = rotations_of(cell);

  auto const mesh32 = *kpoint::Mesh::of({32, 32, 32}); // 32768 points
  auto const mesh64 = *kpoint::Mesh::of({64, 64, 64}); // 262144 points

  BENCHMARK("32^3 mesh") {
    return kpoint::ReciprocalMesh::from_rotations(mesh32, rotations,
                                                  TimeReversal::on)
        .num_irreducible();
  };
  BENCHMARK("64^3 mesh") {
    return kpoint::ReciprocalMesh::from_rotations(mesh64, rotations,
                                                  TimeReversal::on)
        .num_irreducible();
  };
}

TEST_CASE("Brillouin-zone relocation", "[!benchmark]") {
  std::mt19937 rng(17);
  Cell const cell = supercell(2, 0.0, rng);
  std::vector<Matrix3i> const rotations = rotations_of(cell);
  auto const mesh = *kpoint::Mesh::of({24, 24, 24});
  auto const reduced =
      kpoint::ReciprocalMesh::from_rotations(mesh, rotations, TimeReversal::on);
  Lattice const reciprocal{Matrix3d::Identity()};

  BENCHMARK("24^3 mesh into the first BZ") {
    return reduced.brillouin_zone(reciprocal).addresses().size();
  };
}

TEST_CASE("random structure generation", "[!benchmark]") {
  // Fm-3m (hall 523) and P2_1/c (hall 81): a high-symmetry cubic setting where
  // orbits are large and few placements are needed, and a common low-symmetry
  // one where the attempt loop runs longer.
  auto const &cubic = group::SpaceGroup::of(*HallNumber::of(GroupFamily::space, 523));
  auto const &monoclinic = group::SpaceGroup::of(*HallNumber::of(GroupFamily::space, 81));

  generate::Composition const comp{{8, 4}, {14, 8}}; // O4 Si8
  generate::GenerateOptions options;
  options.seed = 2024;

  generate::Generator<group::SpaceGroup> const cubic_gen(cubic, options);
  generate::Generator<group::SpaceGroup> const monoclinic_gen(monoclinic,
                                                              options);

  BENCHMARK("Fm-3m, O4Si8") { return cubic_gen(comp).has_value(); };
  BENCHMARK("P2_1/c, O4Si8") { return monoclinic_gen(comp).has_value(); };
}
