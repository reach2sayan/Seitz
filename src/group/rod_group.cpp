#include <spglib/group/rod_group.hpp>

#include <spglib/data/rod_database.hpp>
#include <spglib/math/fractional.hpp>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <ranges>
#include <span>
#include <vector>

namespace spglib::group {

namespace {

constexpr double kTol = 1e-6;

// --- Linear-algebra helpers (the affine generalisation of the point-group
// subspace arrangement in point_group.cpp) ----------------------------------

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

[[nodiscard]] Eigen::MatrixXd null_space(Eigen::MatrixXd const &m) {
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeFullV);
  Eigen::VectorXd const sv = svd.singularValues();
  auto const rank =
      std::ranges::count_if(sv, [](double x) { return x > kTol; });
  return svd.matrixV().rightCols(m.cols() - rank);
}

// Orthonormal basis of the intersection of two column spaces.
[[nodiscard]] Eigen::MatrixXd intersect_dirs(Eigen::MatrixXd const &a,
                                             Eigen::MatrixXd const &b) {
  if (a.cols() == 0 || b.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  Eigen::MatrixXd stacked(3, a.cols() + b.cols());
  stacked.leftCols(a.cols()) = a;
  stacked.rightCols(b.cols()) = -b;
  Eigen::MatrixXd const ker = null_space(stacked);
  if (ker.cols() == 0) {
    return Eigen::MatrixXd(3, 0);
  }
  return column_space(a * ker.topRows(a.cols()));
}

// --- Affine loci (a point + a direction subspace), periodic along one axis ---

struct RodLocus {
  Vector3d point{Vector3d::Zero()};
  Eigen::MatrixXd dir{Eigen::MatrixXd(3, 0)};
  [[nodiscard]] int dim() const { return static_cast<int>(dir.cols()); }
};

[[nodiscard]] Vector3d unit_axis(int axis) {
  Vector3d e = Vector3d::Zero();
  e[axis] = 1.0;
  return e;
}

// Fold only the periodic component into [0, 1).
[[nodiscard]] Vector3d reduce_periodic(Vector3d p, int axis) {
  p[axis] = math::wrap_to_unit_cell(p[axis]);
  return p;
}

// Equality modulo the rod lattice: exact on the aperiodic axes, mod 1 on the
// periodic axis.
[[nodiscard]] bool same_point_mod(Vector3d const &a, Vector3d const &b,
                                  int axis) {
  Vector3d d = a - b;
  d[axis] -= std::round(d[axis]);
  return d.cwiseAbs().maxCoeff() < kTol;
}

[[nodiscard]] Matrix3d projector(Eigen::MatrixXd const &dir) {
  if (dir.cols() == 0) {
    return Matrix3d::Zero();
  }
  return dir * dir.transpose();
}

// Re-express a direction subspace as an axis-separated basis: the in-plane
// (aperiodic) directions first, then the periodic axis itself if the subspace
// contains it. SVD bases can mix the periodic axis into the in-plane vectors;
// this canonical form lets the generator sample each free coordinate correctly
// (full repeat along the periodic axis, a centred band across the section).
[[nodiscard]] Eigen::MatrixXd canonicalize_dir(Eigen::MatrixXd const &dir,
                                               int axis) {
  Vector3d const e = unit_axis(axis);
  bool const has_periodic =
      dir.cols() > 0 &&
      (dir * (dir.transpose() * e) - e).cwiseAbs().maxCoeff() < kTol;

  // In-plane part: zero the periodic component of every direction, then take a
  // basis of what remains.
  Eigen::MatrixXd projected = dir;
  for (int c = 0; c < projected.cols(); ++c) {
    projected(axis, c) = 0.0;
  }
  Eigen::MatrixXd const plane = column_space(projected);

  int const cols = static_cast<int>(plane.cols()) + (has_periodic ? 1 : 0);
  Eigen::MatrixXd out(3, cols);
  if (plane.cols() > 0) {
    out.leftCols(plane.cols()) = plane;
  }
  if (has_periodic) {
    out.col(cols - 1) = e;
  }
  return out;
}

[[nodiscard]] bool same_locus(RodLocus const &a, RodLocus const &b, int axis) {
  if (a.dim() != b.dim()) {
    return false;
  }
  Matrix3d const pa = projector(a.dir);
  if ((pa - projector(b.dir)).cwiseAbs().maxCoeff() > kTol) {
    return false;
  }
  for (int n = -1; n <= 1; ++n) {
    Vector3d delta = a.point - b.point - static_cast<double>(n) * unit_axis(axis);
    Vector3d const perp = delta - pa * delta;
    if (perp.cwiseAbs().maxCoeff() < kTol) {
      return true;
    }
  }
  return false;
}

void add_unique(std::vector<RodLocus> &loci, RodLocus const &candidate,
                int axis) {
  if (std::ranges::none_of(loci, [&](RodLocus const &l) {
        return same_locus(l, candidate, axis);
      })) {
    loci.push_back(candidate);
  }
}

// The fixed loci of one operation h(p) = R p + t within one period: the affine
// solution set of (R - I) p = -t + n*e_axis over integer shifts n (the periodic
// lattice). A screw / glide with no consistent shift fixes nothing.
[[nodiscard]] std::vector<RodLocus>
fixed_loci(SymmetryOperation const &op, int axis) {
  Matrix3d const m = op.rotation.cast<double>() - Matrix3d::Identity();
  Eigen::MatrixXd const ker = null_space(m); // 3 x k
  Eigen::ColPivHouseholderQR<Matrix3d> const qr(m);

  std::vector<RodLocus> loci;
  for (int n = -1; n <= 2; ++n) {
    Vector3d const rhs = -op.translation + static_cast<double>(n) * unit_axis(axis);
    Vector3d const particular = qr.solve(rhs);
    if ((m * particular - rhs).cwiseAbs().maxCoeff() < kTol) {
      add_unique(loci, RodLocus{reduce_periodic(particular, axis), ker}, axis);
    }
  }
  return loci;
}

// All affine intersections of two loci, trying the neighbouring periodic shifts
// (so a c-free line can meet a c-fixed plane at every period offset).
[[nodiscard]] std::vector<RodLocus>
intersect(RodLocus const &a, RodLocus const &b, int axis) {
  std::vector<RodLocus> results;
  int const ka = a.dim();
  int const kb = b.dim();
  for (int shift = -1; shift <= 1; ++shift) {
    Vector3d const p2 = b.point + static_cast<double>(shift) * unit_axis(axis);
    Vector3d const rhs = p2 - a.point;
    if (ka + kb == 0) {
      if (rhs.cwiseAbs().maxCoeff() < kTol) {
        results.push_back(RodLocus{reduce_periodic(a.point, axis),
                                   Eigen::MatrixXd(3, 0)});
      }
      continue;
    }
    Eigen::MatrixXd stacked(3, ka + kb);
    if (ka > 0) {
      stacked.leftCols(ka) = a.dir;
    }
    if (kb > 0) {
      stacked.rightCols(kb) = -b.dir;
    }
    Eigen::VectorXd const w = stacked.colPivHouseholderQr().solve(rhs);
    if ((stacked * w - rhs).cwiseAbs().maxCoeff() < kTol) {
      Vector3d common = a.point;
      if (ka > 0) {
        common += a.dir * w.head(ka);
      }
      results.push_back(RodLocus{reduce_periodic(common, axis),
                                 intersect_dirs(a.dir, b.dir)});
    }
  }
  return results;
}

[[nodiscard]] Vector3d generic_point(RodLocus const &l) {
  static constexpr std::array<double, 3> kSeed = {0.1357, 0.2468, 0.3791};
  auto const idx = std::views::iota(0, l.dim());
  return std::accumulate(idx.begin(), idx.end(), Vector3d{l.point},
                         [&](Vector3d const &acc, int i) -> Vector3d {
                           return acc + kSeed[static_cast<std::size_t>(i)] *
                                            l.dir.col(i);
                         });
}

// The image of a locus under an operation (a locus of the same dimension; the
// direction basis is re-orthonormalised because the rotation need not be
// orthogonal in the conventional basis).
[[nodiscard]] RodLocus apply_op(SymmetryOperation const &op,
                                RodLocus const &l, int axis) {
  Eigen::MatrixXd const ndir =
      l.dim() == 0 ? Eigen::MatrixXd(3, 0)
                   : column_space(op.rotation.cast<double>() * l.dir);
  return RodLocus{reduce_periodic(op.apply(l.point), axis), ndir};
}

// One derived Wyckoff position before letter assignment.
struct Derived {
  int multiplicity;
  int dof;
  Vector3d origin;
  Matrix3d locus_basis;
  SymmetryOperations orbit_ops;
  SymmetryOperations site_symmetry;
};

} // namespace

Vector3d RodWyckoff::sample(std::span<double const> params) const {
  auto const idx =
      std::views::iota(0, std::min(dof_, static_cast<int>(params.size())));
  return std::accumulate(idx.begin(), idx.end(), Vector3d{locus_origin_},
                         [&](Vector3d const &acc, int i) -> Vector3d {
                           return acc + params[static_cast<std::size_t>(i)] *
                                            locus_basis_.col(i);
                         });
}

Result<RodGroup> RodGroup::from_number(int number) {
  if (number < 1 || number > data::num_rod_groups()) {
    return leaf::new_error(e_message{
        "RodGroup::from_number: rod-group number out of range "
        "(expected 1..75)"});
  }

  SymmetryOperations ops = data::rod_operations_from_database(number);
  if (ops.empty()) {
    return leaf::new_error(e_message{
        "RodGroup::from_number: no operations for this rod group (is the "
        "generated rod table present? run tools/transcribe_rod_groups.py)"});
  }

  RodGroup rg;
  rg.number_ = number;
  rg.symbol_ = data::rod_symbol(number);
  rg.periodic_axis_ = data::rod_periodic_axis();
  rg.operations_ = std::move(ops);
  int const axis = rg.periodic_axis_;

  // 1. The arrangement of candidate loci: the whole cell (general position),
  //    every operation's fixed loci, and the closure under affine intersection.
  std::vector<RodLocus> loci;
  add_unique(loci, RodLocus{Vector3d::Zero(), Matrix3d::Identity()}, axis);
  for (auto const &op : rg.operations_) {
    for (auto const &l : fixed_loci(op, axis)) {
      add_unique(loci, l, axis);
    }
  }
  for (std::size_t i = 0; i < loci.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) {
      for (auto const &l : intersect(loci[i], loci[j], axis)) {
        add_unique(loci, l, axis);
      }
    }
  }

  // 2. One Wyckoff position per orbit of loci under the group (conjugate
  //    stabilisers). For the representative locus, compute its orbit (distinct
  //    images of a generic point, mod lattice), site symmetry and coset
  //    representatives — as SpaceGroup::build_position, with the periodic axis
  //    folded and the aperiodic axes left raw.
  std::vector<bool> assigned(loci.size(), false);
  std::vector<Derived> derived;
  for (std::size_t i = 0; i < loci.size(); ++i) {
    if (assigned[i]) {
      continue;
    }
    RodLocus const &locus = loci[i];
    Vector3d const p0 = generic_point(locus);
    SymmetryOperations site_symmetry;
    SymmetryOperations orbit_ops;
    std::vector<Vector3d> orbit_points;
    for (auto const &op : rg.operations_) {
      Vector3d const image = reduce_periodic(op.apply(p0), axis);
      if (same_point_mod(image, p0, axis)) {
        site_symmetry.push_back(op);
      }
      if (std::ranges::none_of(orbit_points, [&](Vector3d const &q) {
            return same_point_mod(q, image, axis);
          })) {
        orbit_points.push_back(image);
        orbit_ops.push_back(op);
      }
    }
    // Zero-pad the axis-separated direction basis into a fixed Matrix3d (only
    // the first dof columns are meaningful), matching PointWyckoff's storage.
    Eigen::MatrixXd const cdir = canonicalize_dir(locus.dir, axis);
    Matrix3d locus_basis = Matrix3d::Zero();
    if (cdir.cols() > 0) {
      locus_basis.leftCols(cdir.cols()) = cdir;
    }
    derived.push_back(Derived{static_cast<int>(orbit_points.size()),
                              locus.dim(), reduce_periodic(locus.point, axis),
                              locus_basis, std::move(orbit_ops),
                              std::move(site_symmetry)});

    // Mark every locus in this one's orbit as the same Wyckoff position.
    for (auto const &op : rg.operations_) {
      RodLocus const mapped = apply_op(op, locus, axis);
      for (std::size_t k = 0; k < loci.size(); ++k) {
        if (!assigned[k] && same_locus(loci[k], mapped, axis)) {
          assigned[k] = true;
        }
      }
    }
  }

  // 3. Assign letters ('a' = most special; the general position ends up last).
  std::ranges::sort(derived, [](Derived const &x, Derived const &y) {
    if (x.multiplicity != y.multiplicity) {
      return x.multiplicity < y.multiplicity;
    }
    return x.dof < y.dof;
  });
  rg.positions_.reserve(derived.size());
  char letter = 'a';
  for (auto &d : derived) {
    rg.positions_.push_back(RodWyckoff{d.multiplicity, d.dof, letter,
                                       std::move(d.origin),
                                       std::move(d.locus_basis),
                                       std::move(d.orbit_ops),
                                       std::move(d.site_symmetry)});
    ++letter;
  }

  return rg;
}

} // namespace spglib::group
