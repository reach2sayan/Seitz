#include <cppcrystal/spin/spin.hpp>

#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/math/fractional.hpp> // math::nearest_offset
#include <cppcrystal/symmetry/find_symmetry.hpp> // is_overlap_same_type
#include <cppcrystal/symmetry/primitive.hpp>

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
                                        bool time_reversal,
                                        bool with_time_reversal,
                                        bool is_axial) {
    double dst = (with_time_reversal && time_reversal) ? -src : src;
    if (is_axial) {
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
                                          bool time_reversal,
                                          bool with_time_reversal,
                                          bool is_axial) {
    Vector3d dst = rot_cart * v;
    if (with_time_reversal && time_reversal) {
      dst = -dst;
    }
    if (is_axial) {
      dst *= rot_cart.determinant();
    }
    return dst;
  }
  [[nodiscard]] static bool close(Vector3d const &a, Vector3d const &b,
                                  double tol) {
    return (a - b).cwiseAbs().maxCoeff() < tol;
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
                            Matrix3d const &rot_cart, bool with_time_reversal,
                            bool is_axial, double mag_symprec) {
  for (int timerev = 0; timerev <= 1; ++timerev) {
    M const transformed = MomentOps<M>::transform(
        m_j, rot_cart, timerev != 0, with_time_reversal, is_axial);
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
  std::ranges::transform(operations, std::back_inserter(out),
                         [&](auto const &op) {
                           return rotation_in_cartesian(lattice, op.rotation);
                         });
  return out;
}

// The core filter: for each spatial operation, check that it consistently maps
// every site's moment onto the moment of the site it overlaps, determining the
// spin-flip sign. Undetermined operations (all touched moments zero) are kept
// as ordinary, or — with time reversal — as both an ordinary and an
// anti-operation.
template <class M>
[[nodiscard]] MagneticSymmetryOperations
get_operations(SymmetryOperations const &sym_nonspin, MagneticCell const &mcell,
               bool with_time_reversal, bool is_axial, double symprec,
               double mag_symprec) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const rot_cart = cartesian_rotations(cell.lattice(), sym_nonspin);

  MagneticSymmetryOperations out;
  out.reserve(2 * sym_nonspin.size()); // upper bound: 1–2 ops emitted per spatial op
  for (auto const &[rc, op] : std::views::zip(rot_cart, sym_nonspin)) {
    bool found = true;
    bool determined = false;
    int sign = 0;

    for (Index j = 0; j < n; ++j) {
      Vector3d const pos = op.apply(cell.position(j));
      auto const image = std::views::iota(Index{0}, n);
      auto const k_it = std::ranges::find_if(image, [&](Index k) {
        return is_overlap_same_type(cell.position(k), pos, cell.type(k),
                                    cell.type(j), cell.lattice(), symprec);
      });
      if (k_it == image.end()) {
        // Rare: failure to overlap (e.g. too loose symprec); skip the op.
        found = false;
        break;
      }
      Index const k = *k_it;

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

      int const s =
          spin_sign(m_j, m_k, rc, with_time_reversal, is_axial, mag_symprec);

      if (!determined) {
        sign = s;
        determined = true;
        if (sign == 0 || (!with_time_reversal && sign != 1)) {
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
      out.push_back({op.rotation, op.translation, sign == -1});
    } else if (with_time_reversal) {
      out.push_back({op.rotation, op.translation, false}); // sign = +1
      out.push_back({op.rotation, op.translation, true});  // sign = -1
    } else {
      out.push_back({op.rotation, op.translation, false});
    }
  }
  return out;
}

// Flat permutation buffer: row p, entry i = image of atom i under operation p,
// matching both overlap and transformed moment. nullopt is unreachable in
// theory (every site must map somewhere).
template <class M>
[[nodiscard]] std::optional<std::vector<int>>
get_permutations(MagneticSymmetryOperations const &operations,
                 MagneticCell const &mcell, bool with_time_reversal,
                 bool is_axial, double symprec, double mag_symprec) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const rot_cart = cartesian_rotations(cell.lattice(), operations);

  std::vector<int> perm;
  perm.reserve(operations.size() * static_cast<std::size_t>(n));
  for (auto const &[op, rc] : std::views::zip(operations, rot_cart)) {
    for (Index i = 0; i < n; ++i) {
      Vector3d const pos = op.spatial().apply(cell.position(i));
      M const moment =
          MomentOps<M>::transform(MomentOps<M>::moment(mcell, i), rc,
                                  op.time_reversal, with_time_reversal,
                                  is_axial);

      auto const image = std::views::iota(Index{0}, n);
      auto const j_it = std::ranges::find_if(image, [&](Index j) {
        return is_overlap_same_type(pos, cell.position(j), cell.type(i),
                                    cell.type(j), cell.lattice(), symprec) &&
               MomentOps<M>::close(MomentOps<M>::moment(mcell, j), moment,
                                   mag_symprec);
      });
      if (j_it == image.end()) {
        return std::nullopt;
      }
      perm.push_back(static_cast<int>(*j_it));
    }
  }
  return perm;
}

// equivalent_atoms[i] = representative of i's orbit under the permutations.
[[nodiscard]] std::vector<int>
get_orbits(std::vector<int> const &permutations, std::size_t num_sym,
           Index num_atoms) {
  auto const n = static_cast<std::size_t>(num_atoms);
  std::vector<std::optional<int>> equivalent(n);
  for (std::size_t i = 0; i < n; ++i) {
    if (equivalent[i]) {
      continue;
    }
    equivalent[i] = static_cast<int>(i);
    for (std::size_t s = 0; s < num_sym; ++s) {
      auto const image = static_cast<std::size_t>(permutations[s * n + i]);
      equivalent[image] = static_cast<int>(i);
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
                    MagneticCell const &mcell, bool with_time_reversal,
                    bool is_axial) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  auto const &operations = search.operations;
  auto const rot_cart = cartesian_rotations(cell.lattice(), operations);

  // inv_perm[p][i] = the atom that operation p maps onto atom i.
  std::vector<std::vector<int>> inv_perm;
  inv_perm.reserve(operations.size());
  for (auto const row : search.permutation_rows()) {
    std::vector<int> inv(static_cast<std::size_t>(n));
    for (auto const [j, image] : row | std::views::enumerate) {
      inv[static_cast<std::size_t>(image)] = static_cast<int>(j);
    }
    inv_perm.push_back(std::move(inv));
  }

  Positions positions(n, 3);
  std::vector<M> moments(static_cast<std::size_t>(n));

  auto const denom = static_cast<double>(operations.size());
  for (Index i = 0; i < n; ++i) {
    Vector3d pos_res = Vector3d::Zero();
    M moment_res = MomentOps<M>::zero();
    for (auto const &[op, rc, inv_row] :
         std::views::zip(operations, rot_cart, inv_perm)) {
      Index const j = inv_row[static_cast<std::size_t>(i)];
      Vector3d const pos_tmp = op.spatial().apply(cell.position(j));
      // Subtract the input position so the accumulated residual stays small;
      // the per-component nint removes the lattice translation.
      Vector3d const diff = pos_tmp - cell.position(i);
      pos_res += math::nearest_offset(diff);
      moment_res += MomentOps<M>::transform(MomentOps<M>::moment(mcell, j), rc,
                                            op.time_reversal,
                                            with_time_reversal, is_axial) -
                    MomentOps<M>::moment(mcell, i);
    }
    positions.row(i) = (cell.position(i) + pos_res / denom).transpose();
    moments[static_cast<std::size_t>(i)] =
        MomentOps<M>::moment(mcell, i) + moment_res / denom;
  }

  return MagneticCell(Cell(cell.lattice(), positions, cell.types()),
                      MomentOps<M>::pack(std::move(moments)));
}

} // namespace

std::vector<Vector3d>
collect_pure_translations(MagneticSymmetryOperations const &operations) {
  std::vector<Vector3d> out;
  std::ranges::copy(
      operations | std::views::filter([](auto const &op) {
        return op.rotation == Matrix3i::Identity() && !op.time_reversal;
      }) | std::views::transform(&MagneticSymmetryOperation::translation),
      std::back_inserter(out));
  return out;
}

MagneticCell idealized_cell(MagneticSymmetrySearch const &search,
                            MagneticCell const &mcell, bool with_time_reversal,
                            bool is_axial) {
  return visit_moment_kind(mcell, [&]<class M>() {
    return idealized_cell_impl<M>(search, mcell, with_time_reversal, is_axial);
  });
}

Result<MagneticSymmetrySearch>
operations_with_site_tensors(SymmetryOperations const &sym_nonspin,
                             MagneticCell const &mcell, bool with_time_reversal,
                             bool is_axial, double symprec,
                             AngleTolerance angle_tolerance,
                             std::optional<double> mag_symprec) {
  (void)angle_tolerance; // primitive_lattice_vectors uses the Delaunay path
  double const mag_tol = mag_symprec.value_or(symprec);

  MagneticSymmetryOperations operations =
      visit_moment_kind(mcell, [&]<class M>() {
        return get_operations<M>(sym_nonspin, mcell, with_time_reversal,
                                 is_axial, symprec, mag_tol);
      });

  auto permutations = visit_moment_kind(mcell, [&]<class M>() {
    return get_permutations<M>(operations, mcell, with_time_reversal, is_axial,
                               symprec, mag_tol);
  });
  if (!permutations) {
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }
  std::vector<int> equivalent_atoms =
      get_orbits(*permutations, operations.size(), mcell.size());

  std::vector<Vector3d> const pure_trans = collect_pure_translations(operations);
  auto prim_lattice = symmetry::primitive_lattice_vectors(mcell.cell(),
                                                          pure_trans, symprec);
  if (!prim_lattice) {
    // By definition the number of pure translations must be preserved; failure
    // to span the primitive lattice means it was not.
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }

  return MagneticSymmetrySearch{std::move(operations),
                                std::move(equivalent_atoms),
                                std::move(*permutations), *prim_lattice};
}

} // namespace cppcrystal::spin
