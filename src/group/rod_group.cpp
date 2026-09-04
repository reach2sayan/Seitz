#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/group/rod_group.hpp>

#include "data/rod_database.hpp"
#include "group/locus_arrangement.hpp"
#include "math/fractional.hpp"
#include "math/subspace.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cmath>
#include <ranges>
#include <span>
#include <vector>

namespace cppcrystal::group {

namespace {

constexpr double kTol = 1e-6;

[[nodiscard]] Vector3d unit_axis(int axis) {
  Vector3d e = Vector3d::Zero();
  e[axis] = 1.0;
  return e;
}

// The geometry policy of a rod group: loci are affine subspaces (a point plus
// a direction span), periodic along one axis — that component is folded into
// [0, 1) and compared modulo the rod lattice, the aperiodic components are
// exact. Orbit counts come from the distinct images of a generic point.
struct RodGeometry {
  // An affine locus: point + span(dir columns).
  struct Locus {
    Vector3d point{Vector3d::Zero()};
    Eigen::MatrixXd dir{Eigen::MatrixXd(3, 0)};
    [[nodiscard]] int dim() const { return static_cast<int>(dir.cols()); }
  };

  std::span<SymmetryOperation const> ops;
  int axis = 2;

  // Fold only the periodic component into [0, 1).
  [[nodiscard]] Vector3d reduce(Vector3d p) const {
    p[axis] = math::wrap_to_unit_cell(p[axis]);
    return p;
  }

  // Equality modulo the rod lattice: exact on the aperiodic axes, mod 1 on
  // the periodic axis.
  [[nodiscard]] bool same_point(Vector3d const &a, Vector3d const &b) const {
    Vector3d d = a - b;
    d[axis] -= std::round(d[axis]);
    return approx_zero(d, kTol);
  }

  [[nodiscard]] Locus whole_space() const {
    return Locus{Vector3d::Zero(), Matrix3d::Identity()};
  }

  // The fixed loci of one operation h(p) = R p + t within one period: the
  // affine solution set of (R - I) p = -t + n*e_axis over integer shifts n
  // (the periodic lattice). A screw / glide with no consistent shift fixes
  // nothing. Duplicate shifts are deduplicated by the arrangement driver.
  [[nodiscard]] std::vector<Locus>
  fixed_loci(SymmetryOperation const &op) const {
    Matrix3d const m = op.rotation.cast<double>() - Matrix3d::Identity();
    Eigen::MatrixXd const ker = math::null_space(m, kTol); // 3 x k
    Eigen::FullPivLU<Matrix3d> const lu(m);

    std::vector<Locus> loci;
    for (int n = -1; n <= 2; ++n) {
      Vector3d const rhs =
          -op.translation + static_cast<double>(n) * unit_axis(axis);
      Vector3d const particular = lu.solve(rhs);
      if (approx_equal(m * particular, rhs, kTol)) {
        loci.push_back(Locus{reduce(particular), ker});
      }
    }
    return loci;
  }

  // All affine intersections of two loci, trying the neighbouring periodic
  // shifts (so a c-free line can meet a c-fixed plane at every period offset).
  [[nodiscard]] std::vector<Locus> intersections(Locus const &a,
                                                 Locus const &b) const {
    std::vector<Locus> results;
    int const ka = a.dim();
    int const kb = b.dim();
    for (int shift = -1; shift <= 1; ++shift) {
      Vector3d const p2 =
          b.point + static_cast<double>(shift) * unit_axis(axis);
      Vector3d const rhs = p2 - a.point;
      if (ka + kb == 0) {
        if (approx_zero(rhs, kTol)) {
          results.push_back(Locus{reduce(a.point), Eigen::MatrixXd(3, 0)});
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
      if (approx_equal(stacked * w, rhs, kTol)) {
        Vector3d common = a.point;
        if (ka > 0) {
          common += a.dir * w.head(ka);
        }
        results.push_back(Locus{
            reduce(common), math::intersect_column_spaces(a.dir, b.dir, kTol)});
      }
    }
    return results;
  }

  [[nodiscard]] bool same_locus(Locus const &a, Locus const &b) const {
    if (a.dim() != b.dim()) {
      return false;
    }
    Matrix3d const pa = math::projector(a.dir);
    if (!approx_equal(pa, math::projector(b.dir), kTol)) {
      return false;
    }
    for (int n = -1; n <= 1; ++n) {
      Vector3d const delta =
          a.point - b.point - static_cast<double>(n) * unit_axis(axis);
      Vector3d const perp = delta - pa * delta;
      if (approx_zero(perp, kTol)) {
        return true;
      }
    }
    return false;
  }

  // The image of a locus under an operation (a locus of the same dimension;
  // the direction basis is re-orthonormalised because the rotation need not
  // be orthogonal in the conventional basis).
  [[nodiscard]] Locus image(SymmetryOperation const &op, Locus const &l) const {
    Eigen::MatrixXd const ndir =
        l.dim() == 0
            ? Eigen::MatrixXd(3, 0)
            : math::column_space(op.rotation.cast<double>() * l.dir, kTol);
    return Locus{reduce(op.apply(l.point)), ndir};
  }

  [[nodiscard]] static Vector3d generic_point(Locus const &l) {
    static constexpr std::array<double, 3> kSeed = {0.1357, 0.2468, 0.3791};
    auto const idx = std::views::iota(0, l.dim());
    return std::ranges::fold_left(
        idx, Vector3d{l.point}, [&](Vector3d acc, int i) -> Vector3d {
          return acc + kSeed[static_cast<std::size_t>(i)] * l.dir.col(i);
        });
  }

  // Re-express a direction subspace as an axis-separated basis: the in-plane
  // (aperiodic) directions first, then the periodic axis itself if the
  // subspace contains it. SVD bases can mix the periodic axis into the
  // in-plane vectors; this canonical form lets the generator sample each free
  // coordinate correctly (full repeat along the periodic axis, a centred band
  // across the section).
  [[nodiscard]] Eigen::MatrixXd
  canonicalize_dir(Eigen::MatrixXd const &dir) const {
    Vector3d const e = unit_axis(axis);
    bool const has_periodic =
        dir.cols() > 0 && approx_equal(dir * (dir.transpose() * e), e, kTol);

    // In-plane part: zero the periodic component of every direction, then
    // take a basis of what remains.
    Eigen::MatrixXd projected = dir;
    projected.row(axis).setZero();
    Eigen::MatrixXd const plane = math::column_space(projected, kTol);

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

  // Orbit (distinct images of a generic point, mod lattice), site symmetry
  // and coset representatives, with the periodic axis folded and the
  // aperiodic axes left raw.
  [[nodiscard]] detail::DerivedLocus derive(Locus const &locus) const {
    auto partition = detail::partition_orbit(
        ops, generic_point(locus),
        [this](Vector3d const &p) { return reduce(p); },
        [this](Vector3d const &a, Vector3d const &b) {
          return same_point(a, b);
        });
    // Zero-pad the axis-separated direction basis into a fixed Matrix3d (only
    // the first dof columns are meaningful).
    Eigen::MatrixXd const cdir = canonicalize_dir(locus.dir);
    Matrix3d locus_basis = Matrix3d::Zero();
    if (cdir.cols() > 0) {
      locus_basis.leftCols(cdir.cols()) = cdir;
    }
    return detail::DerivedLocus{static_cast<int>(partition.points.size()),
                                locus.dim(),
                                reduce(locus.point),
                                locus_basis,
                                std::move(partition.orbit_ops),
                                std::move(partition.site_symmetry)};
  }
};

} // namespace

Result<RodGroup> RodGroup::from_number(int number) {
  if (number < 1 || number > data::num_rod_groups()) {
    return leaf::new_error(
        e_message{"RodGroup::from_number: rod-group number out of range "
                  "(expected 1..75)"});
  }

  Operations ops = data::rod_operations_from_database(number);
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

  rg.positions_ = detail::derive_wyckoff_positions(
      rg.operations_, RodGeometry{rg.operations_, rg.periodic_axis_});

  return rg;
}

} // namespace cppcrystal::group
