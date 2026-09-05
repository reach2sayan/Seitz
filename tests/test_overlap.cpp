#include "core/overlap.hpp"

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <random>
#include <vector>

using namespace seitz;

namespace {
// Body-centered cubic: atoms at (0,0,0) and (1/2,1/2,1/2).
Cell bcc(double a) {
  Matrix3d lattice = Matrix3d::Identity() * a;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  return Cell(Lattice{lattice}, pos, {0, 0});
}

Matrix3i inversion() { return -Matrix3i::Identity(); }

// The greedy matching the checker replaced, kept verbatim as the reference:
// originals in order, each taking the lowest-index untaken image of its type
// that overlaps it.
bool reference_total_overlap(Cell const &cell, Vector3d const &trans,
                             Matrix3i const &rot, double symprec) {
  Index const n = cell.size();
  Matrix3d const rotd = rot.cast<double>();
  std::vector<bool> found(static_cast<std::size_t>(n), false);
  for (Index io = 0; io < n; ++io) {
    Index ir = 0;
    for (; ir < n; ++ir) {
      if (found[static_cast<std::size_t>(ir)]) {
        continue;
      }
      Vector3d const image = rotd * cell.position(ir) + trans;
      if (cell.type(io) == cell.type(ir) &&
          coincident(cell.position(io), image, cell.lattice().matrix(), symprec,
                     cell.periodicity())) {
        found[static_cast<std::size_t>(ir)] = true;
        break;
      }
    }
    if (ir == n) {
      return false;
    }
  }
  return true;
}

// A k x k x k supercell of a 2-atom rocksalt-like motif, with the second
// species split across two types so that type filtering matters.
Cell supercell(int k, double noise, std::mt19937 &rng) {
  Index const n = 2 * k * k * k;
  Matrix3d const lattice = Matrix3d::Identity() * (4.0 * k);
  Positions pos(n, 3);
  Types types;
  std::uniform_real_distribution<double> jitter(-noise, noise);
  Index row = 0;
  for (int a = 0; a < k; ++a) {
    for (int b = 0; b < k; ++b) {
      for (int c = 0; c < k; ++c) {
        Vector3d const base(a, b, c);
        pos.row(row++) = ((base + Vector3d::Zero()) / k).transpose();
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

Matrix3i random_rotation(std::mt19937 &rng) {
  // A signed permutation matrix: always a lattice symmetry of a cubic cell.
  std::array<int, 3> perm{0, 1, 2};
  std::ranges::shuffle(perm, rng);
  Matrix3i r = Matrix3i::Zero();
  for (int i = 0; i < 3; ++i) {
    r(i, perm[static_cast<std::size_t>(i)]) = (rng() % 2 == 0) ? 1 : -1;
  }
  return r;
}
} // namespace

TEST_CASE("coincident uses the Cartesian minimal image", "[overlap]") {
  Matrix3d lat = Matrix3d::Identity() * 5.0;
  CHECK(coincident(Vector3d(0.0, 0.0, 0.0), Vector3d(1.0, -2.0, 3.0), lat, 1e-5,
                   all_periodic()));
  CHECK_FALSE(coincident(Vector3d(0.0, 0.0, 0.0), Vector3d(0.1, 0.0, 0.0), lat,
                         1e-5, all_periodic()));
  // 0.1 fractional along a=5 -> 0.5 Angstrom > symprec.
}

TEST_CASE("identity is always a symmetry", "[overlap]") {
  OverlapChecker checker(bcc(3.0), 1e-5);
  CHECK(checker.check_total_overlap(Vector3d::Zero(), Matrix3i::Identity()));
}

TEST_CASE("inversion is a symmetry of bcc", "[overlap]") {
  OverlapChecker checker(bcc(3.0), 1e-5);
  CHECK(checker.check_total_overlap(Vector3d::Zero(), inversion()));
}

TEST_CASE("a bogus translation is rejected", "[overlap]") {
  OverlapChecker checker(bcc(3.0), 1e-5);
  CHECK_FALSE(checker.check_total_overlap(Vector3d(0.3, 0.0, 0.0),
                                          Matrix3i::Identity()));
}

TEST_CASE("the bcc centering translation is a symmetry", "[overlap]") {
  // (1/2,1/2,1/2) maps atom 0 -> atom 1 and atom 1 -> atom 0 (mod lattice).
  OverlapChecker checker(bcc(3.0), 1e-5);
  CHECK(checker.check_total_overlap(Vector3d(0.5, 0.5, 0.5),
                                    Matrix3i::Identity()));
}

TEST_CASE("symmetry breaks when the two atoms differ in type", "[overlap]") {
  Matrix3d lattice = Matrix3d::Identity() * 3.0;
  Positions pos(2, 3);
  pos.row(0) << 0.0, 0.0, 0.0;
  pos.row(1) << 0.5, 0.5, 0.5;
  OverlapChecker checker(Cell(Lattice{lattice}, pos, {0, 1}),
                         1e-5); // CsCl-like, distinct types
  // The centering translation now swaps unlike atoms -> not a symmetry.
  CHECK_FALSE(checker.check_total_overlap(Vector3d(0.5, 0.5, 0.5),
                                          Matrix3i::Identity()));
  // Inversion keeps each atom in place -> still a symmetry.
  CHECK(checker.check_total_overlap(Vector3d::Zero(), -Matrix3i::Identity()));
}

TEST_CASE("the indexed checker agrees with the greedy scan it replaced",
          "[overlap]") {
  std::mt19937 rng(2024);
  for (int k : {2, 4, 6}) {
    for (double noise : {0.0, 1e-6, 1e-3}) {
      for (double symprec : {1e-5, 1e-2, 0.3}) {
        Cell const cell = supercell(k, noise, rng);
        OverlapChecker const checker(cell, symprec);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (int trial = 0; trial < 40; ++trial) {
          Matrix3i const rot = random_rotation(rng);
          // Half the trials are genuine candidate translations (atom 0 onto
          // some atom of its type), half are random.
          Vector3d trans;
          if (trial % 2 == 0) {
            Index const target =
                static_cast<Index>(rng() % static_cast<unsigned>(cell.size()));
            trans =
                cell.position(target) - rot.cast<double>() * cell.position(0);
          } else {
            trans = Vector3d(unit(rng), unit(rng), unit(rng));
          }
          INFO("k=" << k << " noise=" << noise << " symprec=" << symprec);
          CHECK(checker.check_total_overlap(trans, rot) ==
                reference_total_overlap(cell, trans, rot, symprec));
        }
      }
    }
  }
}

TEST_CASE("check_total_overlap on a 432-atom supercell", "[!benchmark]") {
  std::mt19937 rng(7);
  Cell const cell = supercell(6, 1e-4, rng);
  OverlapChecker const checker(cell, 1e-3);
  Vector3d const centering(0.5 / 6, 0.5 / 6, 0.5 / 6);
  BENCHMARK("accepted centering translation") {
    return checker.check_total_overlap(centering, Matrix3i::Identity());
  };
  BENCHMARK("rejected translation") {
    return checker.check_total_overlap(Vector3d(0.13, 0.02, 0.41),
                                       Matrix3i::Identity());
  };
}
