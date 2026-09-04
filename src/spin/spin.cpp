#include "spin/search.hpp"
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/tolerance.hpp>

#include "core/position_index.hpp"
#include "math/fractional.hpp" // math::nearest_offset
#include "symmetry/primitive.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <variant>
#include <vector>

namespace cppcrystal::spin {

namespace {

// rot_cart = lattice . rot . lattice^-1: the rotation expressed in Cartesian
// coordinates.
[[nodiscard]] Matrix3d rotation_in_cartesian(Matrix3d const &lattice,
                                             Matrix3i const &rot) {
  return lattice * rot.cast<double>() * lattice.inverse();
}

// How one kind of magnetic moment (collinear scalar / non-collinear 3-vector)
// reads, transforms, compares, and packs back into SiteTensors. Everything
// moment-kind-specific lives here; the search algorithms are generic over M.
template <class M> struct MomentOps;

template <> struct MomentOps<double> {
  [[nodiscard]] static double moment(MagneticCell const &mcell, Index i) {
    return mcell.scalar(i);
  }
  // Time reversal flips the sign; an axial tensor additionally picks up
  // det(R).
  [[nodiscard]] static double transform(double src, Matrix3d const &rot_cart,
                                        bool time_reversal, TimeReversal mode,
                                        TensorKind kind) {
    double dst = (mode == TimeReversal::on && time_reversal) ? -src : src;
    if (kind == TensorKind::axial) {
      dst *= rot_cart.determinant();
    }
    return dst;
  }
  [[nodiscard]] static bool close(double a, double b, double tol) {
    return std::abs(a - b) < tol;
  }
  [[nodiscard]] static double zero() { return 0.0; }
  [[nodiscard]] static SiteTensors pack(std::vector<double> moments) {
    return CollinearTensors{std::move(moments)};
  }
};

template <> struct MomentOps<Vector3d> {
  [[nodiscard]] static Vector3d moment(MagneticCell const &mcell, Index i) {
    return mcell.vector(i);
  }
  [[nodiscard]] static Vector3d transform(Vector3d const &v,
                                          Matrix3d const &rot_cart,
                                          bool time_reversal, TimeReversal mode,
                                          TensorKind kind) {
    Vector3d dst = rot_cart * v;
    if (mode == TimeReversal::on && time_reversal) {
      dst = -dst;
    }
    if (kind == TensorKind::axial) {
      dst *= rot_cart.determinant();
    }
    return dst;
  }
  [[nodiscard]] static bool close(Vector3d const &a, Vector3d const &b,
                                  double tol) {
    return approx_equal(a, b, tol);
  }
  [[nodiscard]] static Vector3d zero() { return Vector3d::Zero(); }
  [[nodiscard]] static SiteTensors pack(std::vector<Vector3d> const &moments) {
    return NoncollinearTensors{to_positions(moments)};
  }
};

// Run `f.operator()<M>()` for the moment kind matching the cell's active
// tensor alternative — the single collinear/non-collinear dispatch point.
template <class F>
[[nodiscard]] decltype(auto) visit_moment_kind(MagneticCell const &mcell,
                                               F &&f) {
  return std::visit(
      [&]<class T>(T const &) {
        if constexpr (std::same_as<T, CollinearTensors>) {
          return f.template operator()<double>();
        } else {
          return f.template operator()<Vector3d>();
        }
      },
      mcell.tensors());
}

// Spin-flip sign in {-1, 0, 1} such that `sign * R(moment_j) == moment_k`; 0
// when the two moments are not related by the operation. sign = 1 - 2*timerev
// for the matching timerev.
template <class M>
[[nodiscard]] int spin_sign(M const &m_j, M const &m_k,
                            Matrix3d const &rot_cart, TimeReversal mode,
                            TensorKind kind, double mag_symprec) {
  for (int timerev = 0; timerev <= 1; ++timerev) {
    M const transformed =
        MomentOps<M>::transform(m_j, rot_cart, timerev != 0, mode, kind);
    if (MomentOps<M>::close(m_k, transformed, mag_symprec)) {
      return 1 - 2 * timerev;
    }
  }
  return 0;
}

[[nodiscard]] std::vector<Matrix3d>
cartesian_rotations(Matrix3d const &lattice, auto const &operations) {
  std::vector<Matrix3d> out;
  out.reserve(operations.size());
  std::ranges::transform(
      operations, std::back_inserter(out), [&](auto const &op) {
        return rotation_in_cartesian(
            lattice,
            OperationTraits<std::remove_cvref_t<decltype(op)>>::spatial(op)
                .rotation);
      });
  return out;
}

// The core filter: for each spatial operation, check that it consistently maps
// every site's moment onto the moment of the site it overlaps, determining the
// spin-flip sign. Undetermined operations (all touched moments zero) are kept
// as ordinary, or — with time reversal — as both an ordinary and an
// anti-operation.
template <class M>
[[nodiscard]] MagneticOperations
get_operations(Operations const &sym_nonspin, MagneticCell const &mcell,
               TimeReversal mode, double symprec, double mag_symprec) {
  TensorKind const kind = mcell.kind();
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const rot_cart =
      cartesian_rotations(cell.lattice().matrix(), sym_nonspin);
  PositionIndex const index(cell, symprec);

  std::vector<MagneticSymmetryOperation> out;
  out.reserve(
      2 * sym_nonspin.size()); // upper bound: 1-2 ops emitted per spatial op
  for (auto const &[rc, op] : std::views::zip(rot_cart, sym_nonspin)) {
    bool found = true;
    bool determined = false;
    int sign = 0;

    for (Index const j : std::views::iota(Index{0}, n)) {
      auto const image =
          index.first_match(op.apply(cell.position(j)), cell.type(j));
      if (!image) {
        // Rare: failure to overlap (e.g. too loose symprec); skip the op.
        found = false;
        break;
      }
      Index const k = *image;

      // Sites whose relevant moments are all (near) zero say nothing about the
      // magnetic symmetry. m and -m coincide within mag_symprec when
      // |m| < 0.5*mag_symprec, so test the moments against 0.5*mag_symprec.
      double const half = 0.5 * mag_symprec;
      M const m_j = MomentOps<M>::moment(mcell, j);
      M const m_k = MomentOps<M>::moment(mcell, k);
      if (MomentOps<M>::close(m_j, MomentOps<M>::zero(), half) &&
          MomentOps<M>::close(m_k, MomentOps<M>::zero(), half)) {
        continue;
      }

      int const s = spin_sign(m_j, m_k, rc, mode, kind, mag_symprec);

      if (!determined) {
        sign = s;
        determined = true;
        if (sign == 0 || (mode == TimeReversal::off && sign != 1)) {
          found = false;
          break;
        }
      } else if (s != sign) {
        found = false;
        break;
      }
    }

    if (!found) {
      continue;
    }
    if (determined) {
      // sign == -1 only occurs with time reversal -> anti-operation.
      out.push_back({op, sign == -1});
    } else if (mode == TimeReversal::on) {
      out.push_back({op, false}); // sign = +1
      out.push_back({op, true});  // sign = -1
    } else {
      out.push_back({op, false});
    }
  }
  return MagneticOperations{std::move(out)};
}

// Flat permutation buffer: row p, entry i = image of atom i under operation p,
// matching both overlap and transformed moment. nullopt is unreachable in
// theory (every site must map somewhere).
template <class M>
[[nodiscard]] std::optional<std::vector<int>>
get_permutations(MagneticOperations const &operations,
                 MagneticCell const &mcell, TimeReversal mode, double symprec,
                 double mag_symprec) {
  TensorKind const kind = mcell.kind();
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const rot_cart =
      cartesian_rotations(cell.lattice().matrix(), operations);
  PositionIndex const index(cell, symprec);

  std::vector<int> perm;
  perm.reserve(operations.size() * static_cast<std::size_t>(n));
  for (auto const &[op, rc] : std::views::zip(operations, rot_cart)) {
    for (Index const i : std::views::iota(Index{0}, n)) {
      Vector3d const pos = op.spatial.apply(cell.position(i));
      M const moment = MomentOps<M>::transform(
          MomentOps<M>::moment(mcell, i), rc, op.time_reversal, mode, kind);

      auto const j = index.first_match(pos, cell.type(i), [&](int k) {
        return MomentOps<M>::close(MomentOps<M>::moment(mcell, k), moment,
                                   mag_symprec);
      });
      if (!j) {
        return std::nullopt;
      }
      perm.push_back(*j);
    }
  }
  return perm;
}

// equivalent_atoms[i] = representative of i's orbit under the permutations
// (perm[p, i] = image of atom i under operation p).
[[nodiscard]] std::vector<int> get_orbits(md::matrix_view<int const> perm) {
  auto const n = static_cast<std::size_t>(perm.extent(1));
  std::vector<std::optional<int>> equivalent(n);
  for (Index const i :
       std::views::iota(Index{0}, static_cast<Index>(perm.extent(1)))) {
    if (equivalent[static_cast<std::size_t>(i)]) {
      continue;
    }
    equivalent[static_cast<std::size_t>(i)] = static_cast<int>(i);
    for (Index const p :
         std::views::iota(Index{0}, static_cast<Index>(perm.extent(0)))) {
      equivalent[static_cast<std::size_t>(perm[p, i])] = static_cast<int>(i);
    }
  }
  std::vector<int> out;
  out.reserve(n);
  std::ranges::transform(equivalent, std::back_inserter(out),
                         [](std::optional<int> const &e) { return *e; });
  return out;
}

template <class M>
[[nodiscard]] MagneticCell
idealized_cell_impl(MagneticSymmetrySearch const &search,
                    MagneticCell const &mcell, TimeReversal mode) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const &operations = search.operations;
  auto const rot_cart =
      cartesian_rotations(cell.lattice().matrix(), operations);

  // inverse[p, i] = the atom that operation p maps onto atom i.
  auto const perm = search.permutations();
  std::vector<int> inverse_buffer(perm.size());
  md::matrix_view<int> const inverse(inverse_buffer.data(), perm.extent(0),
                                     perm.extent(1));
  for (auto const [p, j] : std::views::cartesian_product(
           std::views::iota(Index{0}, static_cast<Index>(perm.extent(0))),
           std::views::iota(Index{0}, static_cast<Index>(perm.extent(1))))) {
    inverse[p, perm[p, j]] = static_cast<int>(j);
  }

  Positions positions(n, 3);
  std::vector<M> moments(static_cast<std::size_t>(n));

  auto const denom = static_cast<double>(operations.size());
  for (Index i = 0; i < n; ++i) {
    Vector3d pos_res = Vector3d::Zero();
    M moment_res = MomentOps<M>::zero();
    for (auto const &[p, op, rc] :
         std::views::zip(std::views::iota(Index{0}, perm.extent(0)), operations,
                         rot_cart)) {
      Index const j = inverse[p, i];
      Vector3d const pos_tmp = op.spatial.apply(cell.position(j));
      // Subtract the input position so the accumulated residual stays small;
      // the per-component nint removes the lattice translation.
      Vector3d const diff = pos_tmp - cell.position(i);
      pos_res += math::nearest_offset(diff);
      moment_res +=
          MomentOps<M>::transform(MomentOps<M>::moment(mcell, j), rc,
                                  op.time_reversal, mode, mcell.kind()) -
          MomentOps<M>::moment(mcell, i);
    }
    positions.row(i) = (cell.position(i) + pos_res / denom).transpose();
    moments[static_cast<std::size_t>(i)] =
        MomentOps<M>::moment(mcell, i) + moment_res / denom;
  }

  return MagneticCell(Cell(cell.lattice(), positions, cell.types()),
                      MomentOps<M>::pack(std::move(moments)), mcell.kind());
}

} // namespace

template <TimeReversal TR>
MagneticCell SpinSearch::idealized(MagneticSymmetrySearch const &search) const {
  MagneticCell const &mcell = cell_;
  return visit_moment_kind(mcell, [&]<class M>() {
    return idealized_cell_impl<M>(search, mcell, TR);
  });
}

template <TimeReversal TR>
Result<MagneticSymmetrySearch> SpinSearch::operations() const {
  Operations const &sym_nonspin = spatial_;
  MagneticCell const &mcell = cell_;
  MagneticTolerance const &tol = tol_;
  constexpr TimeReversal time_reversal = TR;
  // The angle tolerance plays no part here: primitive_lattice_vectors takes
  // the Delaunay path, which is driven by symprec alone.
  double const symprec = tol.symprec;
  double const mag_tol = tol.moment_or_symprec();

  MagneticOperations operations = visit_moment_kind(mcell, [&]<class M>() {
    return get_operations<M>(sym_nonspin, mcell, time_reversal, symprec,
                             mag_tol);
  });

  auto permutations = visit_moment_kind(mcell, [&]<class M>() {
    return get_permutations<M>(operations, mcell, time_reversal, symprec,
                               mag_tol);
  });
  if (!permutations) {
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }
  std::vector<int> equivalent_atoms = get_orbits(md::matrix_view<int const>(
      permutations->data(), static_cast<Index>(operations.size()),
      mcell.size()));

  std::vector<Vector3d> const pure_trans = operations.pure_translations();
  symmetry::PrimitiveFinder<GroupFamily::space> const finder(
      mcell.cell(), {symprec, std::nullopt});
  auto prim_lattice = finder.lattice_from_pure_translations(pure_trans);
  if (!prim_lattice) {
    // By definition the number of pure translations must be preserved; failure
    // to span the primitive lattice means it was not.
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }

  return MagneticSymmetrySearch{
      std::move(operations), std::move(equivalent_atoms),
      std::move(*permutations), prim_lattice->matrix()};
}

template Result<MagneticSymmetrySearch>
SpinSearch::operations<TimeReversal::on>() const;
template Result<MagneticSymmetrySearch>
SpinSearch::operations<TimeReversal::off>() const;
template MagneticCell
SpinSearch::idealized<TimeReversal::on>(MagneticSymmetrySearch const &) const;
template MagneticCell
SpinSearch::idealized<TimeReversal::off>(MagneticSymmetrySearch const &) const;

} // namespace cppcrystal::spin
