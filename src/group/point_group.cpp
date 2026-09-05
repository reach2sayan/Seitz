#include <seitz/core/operation_set.hpp>
#include <seitz/group/point_group.hpp>

#include "core/matrix_order.hpp"
#include "group/locus_arrangement.hpp"
#include "math/subspace.hpp"
#include "symmetry/pointgroup.hpp"
#include <seitz/data/spg_database.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <ranges>
#include <span>
#include <vector>

namespace seitz::group {

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

// Left-coset representatives of the subgroup `stab` in `ops` (one per coset
// g.stab). Applying each to a point of the locus generates its orbit, so the
// count is the multiplicity (|G| / |stab|) exactly — independent of any sampled
// point. Operations have zero translation, so composition is rotation product.
[[nodiscard]] Operations
coset_representatives(std::span<SymmetryOperation const> ops,
                      Operations const &stab) {
  RotationMultimap<int> const by_rotation =
      index_by_rotation(ops, &SymmetryOperation::rotation);
  std::vector<bool> covered(ops.size(), false);
  std::vector<SymmetryOperation> reps;
  for (auto const [i, g] : ops | std::views::enumerate) {
    if (covered[static_cast<std::size_t>(i)]) {
      continue;
    }
    reps.push_back(g);
    for (auto const &h : stab) {
      auto const [lo, hi] = by_rotation.equal_range(g.rotation * h.rotation);
      for (auto const &entry : std::ranges::subrange(lo, hi)) {
        covered[static_cast<std::size_t>(entry.second)] = true;
      }
    }
  }
  return Operations{std::move(reps)};
}

// The geometry policy of a point group: loci are linear subspaces of R^3 held
// as orthonormal bases (dim 0 is the origin), everything is exact about the
// origin, and orbit counts come from cosets rather than sampled points.
struct PointGeometry {
  // A linear subspace, 3 x dim orthonormal basis. The orthogonal projector
  // onto it is cached rather than recomputed: same_locus below is the equality
  // used by add_unique (O(n^2)) and by orbit marking (O(n^2 |G|)), and it used
  // to build both operands' projectors on every one of those comparisons.
  struct Locus {
    explicit Locus(Eigen::MatrixXd b)
        : basis(std::move(b)), projector(math::projector(basis)) {}

    Eigen::MatrixXd basis;
    Matrix3d projector;
    [[nodiscard]] int dim() const { return static_cast<int>(basis.cols()); }
  };

  std::span<SymmetryOperation const> ops;

  [[nodiscard]] Locus whole_space() const {
    return Locus{Matrix3d::Identity()}; // R^3 (general position)
  }

  // The fixed subspace of an operation: the +1 eigenspace, i.e. ker(R - I).
  [[nodiscard]] std::array<Locus, 1>
  fixed_loci(SymmetryOperation const &op) const {
    Eigen::MatrixXd const m = op.rotation.cast<double>() - Matrix3d::Identity();
    return {Locus{math::null_space(m, kTol)}};
  }

  [[nodiscard]] std::array<Locus, 1> intersections(Locus const &a,
                                                   Locus const &b) const {
    return {Locus{math::intersect_column_spaces(a.basis, b.basis, kTol)}};
  }

  [[nodiscard]] bool same_locus(Locus const &a, Locus const &b) const {
    return a.dim() == b.dim() && approx_equal(a.projector, b.projector, kTol);
  }

  // The image of a subspace under an operation (a new subspace of the same
  // dimension; re-orthonormalised because integer rotations need not be
  // Cartesian-orthogonal in the conventional basis).
  [[nodiscard]] Locus image(SymmetryOperation const &op, Locus const &l) const {
    if (l.dim() == 0) {
      return Locus{Eigen::MatrixXd(3, 0)};
    }
    return Locus{
        math::column_space(op.rotation.cast<double>() * l.basis, kTol)};
  }

  // Whether the subspace is fixed pointwise by `rot` (rot v == v for every v
  // in the subspace; trivially true for the origin).
  [[nodiscard]] static bool fixes(Matrix3i const &rot, Locus const &l) {
    Matrix3d const r = rot.cast<double>();
    // Inclusive (`> kTol` negated, i.e. `<= kTol`), so deliberately not
    // approx_equal, which is strict -- see core/tolerance.hpp.
    return std::ranges::none_of(std::views::iota(0, l.dim()), [&](int i) {
      return (r * l.basis.col(i) - l.basis.col(i)).cwiseAbs().maxCoeff() > kTol;
    });
  }

  // The site-symmetry group (the pointwise stabiliser of the locus, exact) and
  // the coset representatives that build the orbit (count = multiplicity).
  [[nodiscard]] detail::DerivedLocus derive(Locus const &l) const {
    Operations const site_symmetry{
        std::from_range,
        ops | std::views::filter([&](SymmetryOperation const &op) {
          return fixes(op.rotation, l);
        })};
    Operations orbit_ops = coset_representatives(ops, site_symmetry);

    Matrix3d locus_basis = Matrix3d::Zero();
    locus_basis.leftCols(l.dim()) = l.basis;

    return detail::DerivedLocus{static_cast<int>(orbit_ops.size()),
                                l.dim(),
                                Vector3d::Zero(),
                                locus_basis,
                                std::move(orbit_ops),
                                std::move(site_symmetry)};
  }
};

} // namespace

Result<PointGroup> PointGroup::from_number(int number) {
  if (number < 1 || number > 32) {
    return leaf::new_error(
        e_message{"PointGroup::from_number: point-group number out of range "
                  "(expected 1..32)"});
  }

  auto const meta = pointgroup_by_number(number);
  int const rep = kRepresentativeSpacegroup[static_cast<std::size_t>(number)];
  std::span<int const> const halls =
      data::halls_with_number<GroupFamily::space>(rep);
  if (halls.empty()) {
    return leaf::new_error(e_message{
        "PointGroup::from_number: representative space group not found"});
  }
  Operations const &conv = data::operations_from_database(
      *HallNumber::of(GroupFamily::space, halls.front()));

  PointGroup pg;
  pg.number_ = number;
  pg.symbol_ = meta.symbol;
  pg.schoenflies_ = meta.schoenflies;
  pg.representative_spacegroup_ = rep;

  // The point group = the distinct rotation parts (zero translation).
  pg.operations_ = Operations{unique_by_rotation(
      conv | std::views::transform([](SymmetryOperation const &op) {
        return SymmetryOperation{.rotation = op.rotation,
                                 .translation = Vector3d::Zero()};
      }),
      &SymmetryOperation::rotation)};

  pg.positions_ = detail::derive_wyckoff_positions(
      pg.operations_, PointGeometry{pg.operations_});

  return pg;
}

} // namespace seitz::group
