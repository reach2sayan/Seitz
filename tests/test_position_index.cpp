// The one property the spatial index must have: for any lattice, tolerance
// and periodicity, its candidate set is a superset of the brute-force
// coincidence set, its match set equals it, and first_match agrees with the
// first hit of a linear scan. Everything built on the index inherits its
// semantics from this.

#include <cppcrystal/core/position_index.hpp>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <algorithm>
#include <optional>
#include <random>
#include <ranges>
#include <vector>

using namespace cppcrystal;

namespace {

struct Scenario {
  char const *name;
  Matrix3d lattice;
  CellPeriodicity periodicity;
};

Matrix3d lattice_from_rows(double a0, double a1, double a2, double b0, double b1,
                           double b2, double c0, double c1, double c2) {
  Matrix3d l;
  // Columns are the basis vectors.
  l.col(0) << a0, a1, a2;
  l.col(1) << b0, b1, b2;
  l.col(2) << c0, c1, c2;
  return l;
}

std::vector<Scenario> scenarios() {
  constexpr CellPeriodicity layer{AxisKind::periodic, AxisKind::periodic,
                                  AxisKind::aperiodic};
  constexpr CellPeriodicity cluster{AxisKind::aperiodic, AxisKind::aperiodic,
                                    AxisKind::aperiodic};
  return {
      {"cubic", Matrix3d::Identity() * 5.0, all_periodic()},
      {"triclinic",
       lattice_from_rows(4.0, 0.3, 0.1, -1.2, 5.5, 0.4, 0.7, -0.9, 6.3),
       all_periodic()},
      {"very skewed",
       lattice_from_rows(1.0, 0.0, 0.0, 0.97, 0.05, 0.0, 0.9, 0.9, 0.02),
       all_periodic()},
      {"layer", lattice_from_rows(3.0, 0.0, 0.0, 1.5, 2.6, 0.0, 0.0, 0.0, 20.0),
       layer},
      {"cluster", Matrix3d::Identity(), cluster},
  };
}

// n random points; every fourth one a near-copy of an earlier point so that
// coincidences at every tolerance actually occur.
Positions random_points(std::mt19937 &rng, Index n, CellPeriodicity const &p,
                        double symprec) {
  std::uniform_real_distribution<double> unit(0.0, 1.0);
  std::uniform_real_distribution<double> wide(-3.0, 3.0);
  std::uniform_real_distribution<double> jitter(-symprec, symprec);
  Positions pos(n, 3);
  for (Index i = 0; i < n; ++i) {
    if (i % 4 == 3) {
      Index const src = static_cast<Index>(rng() % static_cast<unsigned>(i));
      pos.row(i) = pos.row(src);
      for (Index k = 0; k < 3; ++k) {
        pos(i, k) += jitter(rng);
      }
      continue;
    }
    for (Index k = 0; k < 3; ++k) {
      pos(i, k) = p[static_cast<std::size_t>(k)] == AxisKind::periodic
                      ? unit(rng)
                      : wide(rng);
    }
  }
  return pos;
}

std::vector<int> sorted(auto &&range) {
  std::vector<int> out;
  for (int x : range) {
    out.push_back(x);
  }
  std::ranges::sort(out);
  return out;
}

} // namespace

TEST_CASE("PositionIndex agrees with the brute-force scan", "[position_index]") {
  auto const symprec = GENERATE(1e-8, 1e-5, 1e-2, 0.5, 10.0);
  for (auto const &s : scenarios()) {
    INFO(s.name << " symprec=" << symprec);
    std::mt19937 rng(12345);
    constexpr Index n = 200;
    Positions const pos = random_points(rng, n, s.periodicity, symprec);
    Types types(static_cast<std::size_t>(n));
    for (auto &t : types) {
      t = static_cast<int>(rng() % 3);
    }

    PositionIndex const index(
        BucketGeometry::of(s.lattice, symprec, s.periodicity), pos, types,
        s.lattice, symprec, s.periodicity);

    auto const brute = [&](Vector3d const &q) {
      std::vector<int> hits;
      for (Index j = 0; j < n; ++j) {
        if (coincident(q, pos.row(j).transpose(), s.lattice, symprec,
                       s.periodicity)) {
          hits.push_back(static_cast<int>(j));
        }
      }
      return hits;
    };

    // Query at every atom, and at every atom displaced by up to symprec.
    std::uniform_real_distribution<double> jitter(-symprec, symprec);
    for (Index i = 0; i < n; ++i) {
      for (int pass = 0; pass < 2; ++pass) {
        Vector3d q = pos.row(i).transpose();
        if (pass == 1) {
          Vector3d const cart_jitter(jitter(rng), jitter(rng), jitter(rng));
          q += s.lattice.inverse() * cart_jitter;
        }
        auto const expected = brute(q);
        auto const candidates = sorted(index.candidates(q));
        CHECK(std::ranges::includes(candidates, expected));
        CHECK(sorted(index.matches(q)) == expected);

        for (int type = 0; type < 3; ++type) {
          auto typed = expected | std::views::filter([&](int j) {
                         return types[static_cast<std::size_t>(j)] == type;
                       });
          std::optional<int> const first =
              typed.begin() == typed.end() ? std::nullopt
                                           : std::optional<int>(*typed.begin());
          CHECK(index.first_match(q, type) == first);
          CHECK(sorted(index.matches(q, type)) == sorted(typed));
        }
      }
    }
  }
}

TEST_CASE("first_match honours the accept predicate in index order",
          "[position_index]") {
  Matrix3d const lattice = Matrix3d::Identity() * 4.0;
  Positions pos(3, 3);
  pos.row(0) << 0.25, 0.25, 0.25;
  pos.row(1) << 0.25, 0.25, 0.25;
  pos.row(2) << 0.25, 0.25, 0.25;
  Types const types{7, 7, 7};
  PositionIndex const index(BucketGeometry::of(lattice, 1e-5, all_periodic()),
                            pos, types, lattice, 1e-5, all_periodic());
  Vector3d const q(0.25, 0.25, 1.25); // one lattice vector away
  CHECK(index.first_match(q, 7) == 0);
  CHECK(index.first_match(q, 7, [](int j) { return j > 0; }) == 1);
  CHECK(index.first_match(q, 7, [](int j) { return j > 5; }) == std::nullopt);
  CHECK(index.first_match(q, 8) == std::nullopt);
}
