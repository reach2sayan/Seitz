#include <spglib/spin/spin.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/math/fractional.hpp> // math::nint
#include <spglib/symmetry/find_symmetry.hpp> // is_overlap_same_type
#include <spglib/symmetry/primitive.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

namespace spglib::spin {

namespace {

// rot_cart = lattice . rot . lattice^-1: the rotation expressed in Cartesian
// coordinates (spin.c set_rotations_in_cartesian).
[[nodiscard]] Matrix3d rotation_in_cartesian(Matrix3d const &lattice,
                                             Matrix3i const &rot) {
  return lattice * rot.cast<double>() * lattice.inverse();
}

// x -> rot . x + trans (spin.c apply_symmetry_to_position).
[[nodiscard]] Vector3d apply_to_position(Matrix3i const &rot,
                                         Vector3d const &trans,
                                         Vector3d const &x) {
  return rot.cast<double>() * x + trans;
}

// A collinear (scalar) moment transformed by an operation (spin.c
// apply_symmetry_to_site_scalar): time reversal flips the sign; an axial tensor
// additionally picks up det(R).
[[nodiscard]] double apply_to_scalar(double src, Matrix3d const &rot_cart,
                                     bool time_reversal,
                                     bool with_time_reversal, bool is_axial) {
  double dst = (with_time_reversal && time_reversal) ? -src : src;
  if (is_axial) {
    dst *= rot_cart.determinant();
  }
  return dst;
}

// A non-collinear (vector) moment transformed by an operation (spin.c
// apply_symmetry_to_site_vector).
[[nodiscard]] Vector3d apply_to_vector(Vector3d const &v,
                                       Matrix3d const &rot_cart,
                                       bool time_reversal,
                                       bool with_time_reversal, bool is_axial) {
  Vector3d dst = rot_cart * v;
  if (with_time_reversal && time_reversal) {
    dst = -dst;
  }
  if (is_axial) {
    dst *= rot_cart.determinant();
  }
  return dst;
}

// Spin-flip sign in {-1, 0, 1} such that `sign * R(spin_j) == spin_k`; 0 when
// the two moments are not related by the operation (spin.c
// get_operation_sign_on_scalar). sign = 1 - 2*timerev for the matching timerev.
[[nodiscard]] int sign_on_scalar(double spin_j, double spin_k,
                                 Matrix3d const &rot_cart,
                                 bool with_time_reversal, bool is_axial,
                                 double mag_symprec) {
  for (int timerev = 0; timerev <= 1; ++timerev) {
    double const transformed = apply_to_scalar(
        spin_j, rot_cart, timerev != 0, with_time_reversal, is_axial);
    if (std::abs(spin_k - transformed) < mag_symprec) {
      return 1 - 2 * timerev;
    }
  }
  return 0;
}

[[nodiscard]] int sign_on_vector(Vector3d const &v_j, Vector3d const &v_k,
                                 Matrix3d const &rot_cart,
                                 bool with_time_reversal, bool is_axial,
                                 double mag_symprec) {
  for (int timerev = 0; timerev <= 1; ++timerev) {
    Vector3d const transformed = apply_to_vector(
        v_j, rot_cart, timerev != 0, with_time_reversal, is_axial);
    if ((v_k - transformed).cwiseAbs().maxCoeff() < mag_symprec) {
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

// The core filter (spin.c get_operations): for each spatial operation, check
// that it consistently maps every site's moment onto the moment of the site it
// overlaps, determining the spin-flip sign. Undetermined operations (all
// touched moments zero) are kept as ordinary, or — with time reversal — as both
// an ordinary and an anti-operation.
[[nodiscard]] MagneticSymmetryOperations
get_operations(SymmetryOperations const &sym_nonspin, MagneticCell const &mcell,
               bool with_time_reversal, bool is_axial, double symprec,
               double mag_symprec) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  bool const collinear = mcell.rank() == SiteTensor::collinear;
  auto const rot_cart = cartesian_rotations(cell.lattice(), sym_nonspin);

  MagneticSymmetryOperations out;
  out.reserve(2 * sym_nonspin.size()); // upper bound: 1–2 ops emitted per spatial op
  for (std::size_t i = 0; i < sym_nonspin.size(); ++i) {
    auto const &op = sym_nonspin[i];
    bool found = true;
    bool determined = false;
    int sign = 0;

    for (Index j = 0; j < n; ++j) {
      Vector3d const pos =
          apply_to_position(op.rotation, op.translation, cell.position(j));
      Index k = 0;
      for (; k < n; ++k) {
        if (is_overlap_same_type(cell.position(k), pos, cell.type(k),
                                 cell.type(j), cell.lattice(), symprec)) {
          break;
        }
      }
      if (k == n) {
        // Rare: failure to overlap (e.g. too loose symprec); skip the op.
        found = false;
        break;
      }

      // Sites whose relevant moments are all (near) zero say nothing about the
      // magnetic symmetry. m and -m coincide within mag_symprec when
      // |m| < 0.5*mag_symprec, so test the moments against 0.5*mag_symprec.
      double const half = 0.5 * mag_symprec;
      if (collinear) {
        if (std::abs(mcell.scalar(j)) < half &&
            std::abs(mcell.scalar(k)) < half) {
          continue;
        }
      } else {
        if (mcell.vector(j).cwiseAbs().maxCoeff() < half &&
            mcell.vector(k).cwiseAbs().maxCoeff() < half) {
          continue;
        }
      }

      int const s =
          collinear
              ? sign_on_scalar(mcell.scalar(j), mcell.scalar(k), rot_cart[i],
                               with_time_reversal, is_axial, mag_symprec)
              : sign_on_vector(mcell.vector(j), mcell.vector(k), rot_cart[i],
                               with_time_reversal, is_axial, mag_symprec);

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

// permutations[p * size + i] = image of atom i under operation p, matching both
// overlap and transformed moment (spin.c get_symmetry_permutations). nullopt is
// unreachable in theory (every site must map somewhere).
[[nodiscard]] std::optional<std::vector<int>>
get_permutations(MagneticSymmetryOperations const &operations,
                 MagneticCell const &mcell, bool with_time_reversal,
                 bool is_axial, double symprec, double mag_symprec) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  bool const collinear = mcell.rank() == SiteTensor::collinear;
  auto const rot_cart = cartesian_rotations(cell.lattice(), operations);

  std::vector<int> perm(operations.size() * static_cast<std::size_t>(n), -1);
  for (auto const [pi, pair] :
       std::views::zip(operations, rot_cart) | std::views::enumerate) {
    auto const &[op, rc] = pair;
    auto const p = static_cast<std::size_t>(pi);
    for (Index i = 0; i < n; ++i) {
      Vector3d const pos =
          apply_to_position(op.rotation, op.translation, cell.position(i));
      double scalar = 0.0;
      Vector3d vector = Vector3d::Zero();
      if (collinear) {
        scalar = apply_to_scalar(mcell.scalar(i), rc, op.time_reversal,
                                 with_time_reversal, is_axial);
      } else {
        vector = apply_to_vector(mcell.vector(i), rc, op.time_reversal,
                                 with_time_reversal, is_axial);
      }

      for (Index j = 0; j < n; ++j) {
        if (!is_overlap_same_type(pos, cell.position(j), cell.type(i),
                                  cell.type(j), cell.lattice(), symprec)) {
          continue;
        }
        if (collinear && std::abs(mcell.scalar(j) - scalar) >= mag_symprec) {
          continue;
        }
        if (!collinear &&
            (mcell.vector(j) - vector).cwiseAbs().maxCoeff() >= mag_symprec) {
          continue;
        }
        perm[p * static_cast<std::size_t>(n) + static_cast<std::size_t>(i)] =
            static_cast<int>(j);
        break;
      }
      if (perm[p * static_cast<std::size_t>(n) + static_cast<std::size_t>(i)] ==
          -1) {
        return std::nullopt;
      }
    }
  }
  return perm;
}

// equivalent_atoms[i] = representative of i's orbit under the permutations
// (spin.c get_orbits).
[[nodiscard]] std::vector<int>
get_orbits(std::vector<int> const &permutations, std::size_t num_sym,
           Index num_atoms) {
  std::vector<int> equivalent(static_cast<std::size_t>(num_atoms), -1);
  for (Index i = 0; i < num_atoms; ++i) {
    auto const ui = static_cast<std::size_t>(i);
    if (equivalent[ui] != -1) {
      continue;
    }
    equivalent[ui] = static_cast<int>(i);
    for (std::size_t s = 0; s < num_sym; ++s) {
      auto const image = static_cast<std::size_t>(
          permutations[s * static_cast<std::size_t>(num_atoms) + ui]);
      equivalent[image] = static_cast<int>(i);
    }
  }
  return equivalent;
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

MagneticCell idealized_cell(std::vector<int> const &permutations,
                            MagneticCell const &mcell,
                            MagneticSymmetryOperations const &operations,
                            bool with_time_reversal, bool is_axial) {
  Cell const &cell = mcell.cell();
  Index const n = cell.size();
  std::size_t const num_ops = operations.size();
  bool const collinear = mcell.rank() == SiteTensor::collinear;
  auto const rot_cart = cartesian_rotations(cell.lattice(), operations);

  // inv_perm[p][i] = the atom that operation p maps onto atom i.
  std::vector<std::vector<int>> inv_perm(
      num_ops, std::vector<int>(static_cast<std::size_t>(n)));
  for (std::size_t p = 0; p < num_ops; ++p) {
    for (Index j = 0; j < n; ++j) {
      int const image = permutations[p * static_cast<std::size_t>(n) +
                                     static_cast<std::size_t>(j)];
      inv_perm[p][static_cast<std::size_t>(image)] = static_cast<int>(j);
    }
  }

  Positions positions(n, 3);
  CollinearTensors scalars;
  NoncollinearTensors vectors;
  if (collinear) {
    scalars.resize(static_cast<std::size_t>(n));
  } else {
    vectors.resize(n, 3);
  }

  auto const denom = static_cast<double>(num_ops);
  for (Index i = 0; i < n; ++i) {
    Vector3d pos_res = Vector3d::Zero();
    double scalar_res = 0.0;
    Vector3d vector_res = Vector3d::Zero();
    for (auto const &[op, rc, inv_row] :
         std::views::zip(operations, rot_cart, inv_perm)) {
      Index const j = inv_row[static_cast<std::size_t>(i)];
      Vector3d const pos_tmp =
          apply_to_position(op.rotation, op.translation, cell.position(j));
      // Subtract the input position so the accumulated residual stays small;
      // the per-component nint removes the lattice translation.
      Vector3d const diff = pos_tmp - cell.position(i);
      pos_res += math::nearest_offset(diff);
      if (collinear) {
        double const t = apply_to_scalar(mcell.scalar(j), rc, op.time_reversal,
                                         with_time_reversal, is_axial);
        scalar_res += t - mcell.scalar(i);
      } else {
        Vector3d const t = apply_to_vector(mcell.vector(j), rc, op.time_reversal,
                                           with_time_reversal, is_axial);
        vector_res += t - mcell.vector(i);
      }
    }
    positions.row(i) = (cell.position(i) + pos_res / denom).transpose();
    auto const ui = static_cast<std::size_t>(i);
    if (collinear) {
      scalars[ui] = mcell.scalar(i) + scalar_res / denom;
    } else {
      vectors.row(i) = (mcell.vector(i) + vector_res / denom).transpose();
    }
  }

  SiteTensors tensors = collinear
                            ? SiteTensors{std::move(scalars)}
                            : SiteTensors{std::move(vectors)};
  return MagneticCell(Cell(cell.lattice(), positions, cell.types()),
                      std::move(tensors));
}

Result<MagneticSymmetrySearch>
operations_with_site_tensors(SymmetryOperations const &sym_nonspin,
                             MagneticCell const &mcell, bool with_time_reversal,
                             bool is_axial, double symprec,
                             AngleTolerance angle_tolerance,
                             std::optional<double> mag_symprec) {
  (void)angle_tolerance; // primitive_lattice_vectors uses the Delaunay path
  double const mag_tol = mag_symprec.value_or(symprec);

  MagneticSymmetryOperations operations = get_operations(
      sym_nonspin, mcell, with_time_reversal, is_axial, symprec, mag_tol);

  auto permutations = get_permutations(operations, mcell, with_time_reversal,
                                       is_axial, symprec, mag_tol);
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
    // to span the primitive lattice means it was not (spin.c `multi !=
    // pure_trans->size`).
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }

  return MagneticSymmetrySearch{std::move(operations),
                                std::move(equivalent_atoms),
                                std::move(*permutations), *prim_lattice};
}

} // namespace spglib::spin
