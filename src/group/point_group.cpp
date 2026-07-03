#include <cppcrystal/group/point_group.hpp>

#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/symmetry/pointgroup.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <numeric>
#include <ranges>
#include <span>
#include <vector>

namespace cppcrystal::group {

namespace {

constexpr double kTol = 1e-7;

// Representative symmorphic space group for each point-group number (1..32),
// indexed directly (entry 0 unused). Its conventional rotation set IS the point
// group; symmorphic P groups are chosen so the rotation count equals the point
// group order exactly.
constexpr std::array<int, 33> kRepresentativeSpacegroup = {
    0,   1,   2,   3,   6,   10,  16,  25,  47,  75,  81,
    83,  89,  99,  111, 123, 143, 147, 149, 156, 164, 168,
    174, 175, 177, 183, 187, 191, 195, 200, 207, 215, 221};

// A linear subspace of R^3, held as an orthonormal basis (3 x dim). dim() == 0
// is the origin {0}.
struct Subspace {
  Eigen::MatrixXd basis; // 3 x dim, orthonormal columns
  [[nodiscard]] int dim() const { return static_cast<int>(basis.cols()); }
};

// Orthonormal basis of the column space of `m` (columns with non-negligible
// singular value).
[[nodiscard]] Eigen::MatrixXd column_space(Eigen::MatrixXd const &m) {
  if (m.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU);
  Eigen::VectorXd const sv = svd.singularValues();
  auto const rank =
      std::ranges::count_if(sv, [](double x) { return x > kTol; });
  return svd.matrixU().leftCols(rank);
}

// Orthonormal basis of the null space of `m` (the domain-space directions that
// `m` annihilates).
[[nodiscard]] Eigen::MatrixXd null_space(Eigen::MatrixXd const &m) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullV);
  Eigen::VectorXd const sv = svd.singularValues();
  auto const n = m.cols();
  auto const rank =
      std::ranges::count_if(sv, [](double x) { return x > kTol; });
  auto const nullity = n - rank;
  // The last `nullity` columns of V correspond to the smallest singular values.
  return svd.matrixV().rightCols(nullity);
}

// The fixed subspace of an operation: the +1 eigenspace, i.e. ker(R - I).
[[nodiscard]] Subspace fixed_space(Matrix3i const &rot) {
  Eigen::MatrixXd const m =
      rot.cast<double>() - Matrix3d::Identity();
  return Subspace{null_space(m)};
}

// Intersection of two subspaces: vectors expressible in both column spaces.
[[nodiscard]] Subspace intersect(Subspace const &a, Subspace const &b) {
  if (a.dim() == 0 || b.dim() == 0) {
    return Subspace{Eigen::MatrixXd(3, 0)};
  }
  // v = a.basis * x = b.basis * y  <=>  [a.basis  -b.basis] [x; y] = 0.
  Eigen::MatrixXd stacked(3, a.dim() + b.dim());
  stacked.leftCols(a.dim()) = a.basis;
  stacked.rightCols(b.dim()) = -b.basis;
  Eigen::MatrixXd const ker = null_space(stacked); // (da+db) x k
  if (ker.cols() == 0) {
    return Subspace{Eigen::MatrixXd(3, 0)};
  }
  Eigen::MatrixXd const vectors = a.basis * ker.topRows(a.dim());
  return Subspace{column_space(vectors)};
}

[[nodiscard]] Matrix3d projector(Subspace const &s) {
  if (s.dim() == 0) {
    return Matrix3d::Zero();
  }
  return s.basis * s.basis.transpose();
}

[[nodiscard]] bool same_subspace(Subspace const &a, Subspace const &b) {
  return a.dim() == b.dim() &&
         (projector(a) - projector(b)).cwiseAbs().maxCoeff() < kTol;
}

// Whether the subspace is fixed pointwise by `rot` (rot v == v for every v in
// the subspace; trivially true for the origin).
[[nodiscard]] bool fixes(Matrix3i const &rot, Subspace const &s) {
  Matrix3d const r = rot.cast<double>();
  return std::ranges::none_of(std::views::iota(0, s.dim()), [&](int i) {
    return (r * s.basis.col(i) - s.basis.col(i)).cwiseAbs().maxCoeff() > kTol;
  });
}

// The image of a subspace under an operation (a new subspace of the same
// dimension; re-orthonormalised because integer rotations need not be Cartesian-
// orthogonal in the conventional basis).
[[nodiscard]] Subspace apply(Matrix3i const &rot, Subspace const &s) {
  if (s.dim() == 0) {
    return Subspace{Eigen::MatrixXd(3, 0)};
  }
  Eigen::MatrixXd const mapped = rot.cast<double>() * s.basis;
  return Subspace{column_space(mapped)};
}

// Add `s` to `subs` unless an equal subspace is already present.
void add_unique(std::vector<Subspace> &subs, Subspace const &s) {
  if (std::ranges::none_of(
          subs, [&](Subspace const &t) { return same_subspace(t, s); })) {
    subs.push_back(s);
  }
}

// The full arrangement of stabiliser subspaces: the whole space, every
// operation's fixed space, and the closure of those under pairwise
// intersection.
[[nodiscard]] std::vector<Subspace>
stabilizer_subspaces(std::span<SymmetryOperation const> ops) {
  std::vector<Subspace> subs;
  add_unique(subs, Subspace{Matrix3d::Identity()}); // R^3 (general position)
  for (auto const &op : ops) {
    add_unique(subs, fixed_space(op.rotation));
  }
  // Intersection closure (re-scan including newly added subspaces).
  for (std::size_t i = 0; i < subs.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      add_unique(subs, intersect(subs[i], subs[j]));
    }
  }
  return subs;
}

// Left-coset representatives of the subgroup `stab` in `ops` (one per coset
// g.stab). Applying each to a point of the locus generates its orbit, so the
// count is the multiplicity (|G| / |stab|) exactly — independent of any sampled
// point. Operations have zero translation, so composition is rotation product.
[[nodiscard]] SymmetryOperations
coset_representatives(std::span<SymmetryOperation const> ops,
                      SymmetryOperations const &stab) {
  std::vector<bool> covered(ops.size(), false);
  SymmetryOperations reps;
  for (std::size_t i = 0; i < ops.size(); ++i) {
    if (covered[i]) {
      continue;
    }
    reps.push_back(ops[i]);
    for (auto const &h : stab) {
      Matrix3i const product = ops[i].rotation * h.rotation;
      for (std::size_t j = 0; j < ops.size(); ++j) {
        if (ops[j].rotation == product) {
          covered[j] = true;
        }
      }
    }
  }
  return reps;
}

} // namespace

Vector3d PointWyckoff::sample(std::span<double const> params) const {
  auto const idx =
      std::views::iota(0, std::min(dof_, static_cast<int>(params.size())));
  return std::accumulate(idx.begin(), idx.end(), Vector3d{Vector3d::Zero()},
                         [&](Vector3d const &acc, int i) -> Vector3d {
                           return acc + params[static_cast<std::size_t>(i)] *
                                            locus_basis_.col(i);
                         });
}

Result<PointGroup> PointGroup::from_number(int number) {
  if (number < 1 || number > 32) {
    return leaf::new_error(e_message{
        "PointGroup::from_number: point-group number out of range "
        "(expected 1..32)"});
  }

  auto const meta = symmetry::pointgroup_by_number(number);
  int const rep = kRepresentativeSpacegroup[static_cast<std::size_t>(number)];
  std::span<int const> const halls = data::kCatalog.with_number(rep);
  if (halls.empty()) {
    return leaf::new_error(e_message{
        "PointGroup::from_number: representative space group not found"});
  }
  SymmetryOperations const conv = data::operations_from_database(halls.front());

  PointGroup pg;
  pg.number_ = number;
  pg.symbol_ = meta.symbol;
  pg.schoenflies_ = meta.schoenflies;
  pg.representative_spacegroup_ = rep;

  // The point group = the distinct rotation parts (zero translation).
  for (auto const &op : conv) {
    if (std::ranges::none_of(pg.operations_, [&](SymmetryOperation const &o) {
          return o.rotation == op.rotation;
        })) {
      pg.operations_.push_back(SymmetryOperation{op.rotation, Vector3d::Zero()});
    }
  }

  // Derive the Wyckoff positions from the arrangement of stabiliser subspaces:
  // one position per orbit of subspaces under the group.
  std::vector<Subspace> const subs = stabilizer_subspaces(pg.operations_);
  std::vector<bool> assigned(subs.size(), false);

  struct Derived {
    int multiplicity;
    int dof;
    Matrix3d locus_basis;
    SymmetryOperations orbit_ops;
    SymmetryOperations site_symmetry;
  };
  std::vector<Derived> derived;

  for (std::size_t i = 0; i < subs.size(); ++i) {
    if (assigned[i]) {
      continue;
    }
    Subspace const &s = subs[i];

    // The site-symmetry group (the pointwise stabiliser of the locus, exact) and
    // the coset representatives that build the orbit (count = multiplicity).
    SymmetryOperations site_symmetry;
    for (auto const &op : pg.operations_) {
      if (fixes(op.rotation, s)) {
        site_symmetry.push_back(op);
      }
    }
    SymmetryOperations orbit_ops =
        coset_representatives(pg.operations_, site_symmetry);

    Matrix3d locus_basis = Matrix3d::Zero();
    for (int c = 0; c < s.dim(); ++c) {
      locus_basis.col(c) = s.basis.col(c);
    }

    derived.push_back(Derived{static_cast<int>(orbit_ops.size()), s.dim(),
                              locus_basis, std::move(orbit_ops),
                              std::move(site_symmetry)});

    // Mark every subspace in this one's orbit (conjugate stabilisers) as the
    // same Wyckoff position.
    for (auto const &op : pg.operations_) {
      Subspace const mapped = apply(op.rotation, s);
      for (std::size_t k = 0; k < subs.size(); ++k) {
        if (!assigned[k] && same_subspace(subs[k], mapped)) {
          assigned[k] = true;
        }
      }
    }
  }

  // 'a' = the most special (smallest multiplicity); the general position last.
  std::ranges::sort(derived, [](Derived const &x, Derived const &y) {
    if (x.multiplicity != y.multiplicity) {
      return x.multiplicity < y.multiplicity;
    }
    return x.dof < y.dof;
  });

  pg.positions_.reserve(derived.size());
  char letter = 'a';
  for (auto &d : derived) {
    pg.positions_.push_back(PointWyckoff{d.multiplicity, d.dof, letter,
                                         d.locus_basis, std::move(d.orbit_ops),
                                         std::move(d.site_symmetry)});
    ++letter;
  }

  return pg;
}

} // namespace cppcrystal::group
