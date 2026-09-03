#include "magnetic/identify.hpp"

#include "core/matrix_order.hpp"
#include <cppcrystal/core/operation_set.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/msg_database.hpp>
#include "math/fractional.hpp"
#include "math/integer_matrix.hpp"
#include "refine/refinement.hpp"
#include "spacegroup/spacegroup.hpp"
#include "spin/search.hpp"
#include "symmetry/primitive.hpp"

#include <algorithm>
#include <cmath>
#include <ranges>
#include <cstddef>
#include <iterator>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace cppcrystal::magnetic {

namespace {



constexpr int kMaxDenominator = 100;

// Internal construction buffer: operations are appended incrementally, then
// wrapped in a MagneticOperations at the boundary.
using MagOps = std::vector<MagneticSymmetryOperation>;

[[nodiscard]] MagneticSymmetryOperation identity_operation(bool time_reversal) {
  return {{Matrix3i::Identity(), Vector3d::Zero()}, time_reversal};
}

// Whether a magnetic space group is built as the family (ignore time reversal)
// or the maximal (ordinary-only) space group.
enum class SpaceGroupKind { family, maximal };

// Family / maximal space group of a magnetic symmetry, plus its non-magnetic
// symmetry operations (in the input setting).
[[nodiscard]] std::optional<std::pair<Operations, SpacegroupMatch>>
space_group_of_magnetic_symmetry(
    MagneticOperations const &magnetic_symmetry, SpaceGroupKind kind,
    double symprec) {
  bool const ignore_time_reversal = kind == SpaceGroupKind::family;

  // Type-II magnetic space groups contain the pure time-reversal (I, 0)1'.
  bool const is_type2 =
      std::ranges::any_of(magnetic_symmetry, [&](auto const &op) {
        return op.spatial.is_identity_rotation() && op.time_reversal &&
               approx_equal(op.spatial.translation, Vector3d::Zero(), symprec);
      });

  // Keep an operation in the spatial subgroup unless it must be dropped: primed
  // operations are removed from the maximal subgroup, and for a type-II group
  // every g appears as both g and g', so the primed copy is removed.
  auto keep_in_subgroup = [&ignore_time_reversal, &is_type2](const auto &op) {
    bool const drop = (!ignore_time_reversal && op.time_reversal) ||
                      (is_type2 && op.time_reversal);
    return !drop;
  };

  Operations const sym{
      std::from_range,
      magnetic_symmetry | std::views::filter(keep_in_subgroup) |
          std::views::transform([](auto const &op) { return op.spatial; })};

  // (a, b, c) = (a_prim, b_prim, c_prim) . t_mat.
  auto const prim = sym.to_primitive({symprec, std::nullopt});
  if (!prim) {
    return std::nullopt;
  }
  Operations const &prim_sym = prim->first;
  Matrix3d const &t_mat = prim->second;

  // NOT spacegroup_type_from_symmetry: that Niggli-reduces the primitive
  // lattice before matching, which this path must not do, and it discards the
  // t_mat the change of basis below needs. The overlap is superficial.
  auto sg = spacegroup::search_spacegroup_with_symmetry(
      prim_sym, Matrix3d::Identity(), symprec);
  if (!sg) {
    return std::nullopt;
  }
  SpacegroupMatch spacegroup =
      refine::find_similar_bravais_lattice(sg.value(), symprec);
  // Change basis from primitive to original: (a_std,b_std,c_std) = (a,b,c) .
  // t_mat^-1 . P.
  spacegroup.bravais_lattice = t_mat.inverse() * spacegroup.bravais_lattice;
  return std::make_pair(std::move(sym), spacegroup);
}

// Coset representative of XSG in MSG (assumes type III or IV).
// representative[0] is the identity; representative[1] is an anti-translation
// (type IV) or some primed operation (type III).
[[nodiscard]] std::optional<MagOps>
get_representative(MagneticOperations const &magnetic_symmetry) {
  MagOps rep{identity_operation(false)};

  auto const anti_translation =
      std::ranges::find_if(magnetic_symmetry, [](auto const &op) {
        return op.spatial.is_identity_rotation() && op.time_reversal;
      });
  if (anti_translation != magnetic_symmetry.end()) {
    rep.push_back({anti_translation->spatial, true}); // type IV
    return rep;
  }
  auto const primed = std::ranges::find_if(
      magnetic_symmetry, [](auto const &op) { return op.time_reversal; });
  if (primed != magnetic_symmetry.end()) {
    rep.push_back({primed->spatial, true}); // type III
    return rep;
  }
  return std::nullopt;
}

// Index [G:H] of a subgroup from the two group orders, restricted to the values
// that occur between the family/maximal subgroups of a magnetic space group: 1
// (the groups coincide) or 2 (H is a halving subgroup). std::nullopt otherwise.
[[nodiscard]] std::optional<int> subgroup_index(std::size_t order,
                                                std::size_t suborder) {
  if (order == suborder) {
    return 1;
  }
  if (order == 2 * suborder) {
    return 2;
  }
  return std::nullopt;
}

// The MSG type and the coset representatives of XSG in MSG.
//
// The four Opechowski-Guccione types are distinguished by a single index,
// [FSG : XSG]: the family space group FSG (the spatial parts, time reversal
// ignored) over its subgroup XSG of time-even operations.
//   * index 1 - time reversal does not divide the spatial group, so either no
//     operation is primed (type I, colorless) or the group also contains the
//     pure time reversal 1' (type II, grey).
//   * index 2 - XSG is a halving subgroup; the missing spatial operations recur
//     only when primed, as a pure anti-translation (type IV) or a primed point
//     operation (type III), read off the coset representative.
[[nodiscard]] std::optional<std::pair<MagneticType, MagOps>>
magnetic_space_group_type(MagneticOperations const &magnetic_symmetry,
                          std::size_t num_sym_fsg, std::size_t num_sym_xsg) {
  std::size_t const num_sym_msg = magnetic_symmetry.size();
  auto const spatial_index = subgroup_index(num_sym_fsg, num_sym_xsg);
  if (!spatial_index) {
    return std::nullopt;
  }

  if (*spatial_index == 1) {
    if (num_sym_msg == num_sym_fsg) {
      return std::make_pair(MagneticType::type_i,
                            MagOps{identity_operation(false)});
    }
    if (num_sym_msg == 2 * num_sym_fsg) {
      return std::make_pair(
          MagneticType::type_ii,
          MagOps{identity_operation(false), identity_operation(true)});
    }
    return std::nullopt;
  }

  // *spatial_index == 2
  auto rep = get_representative(magnetic_symmetry);
  if (!rep) {
    return std::nullopt;
  }
  MagneticType const type = (*rep)[1].spatial.is_identity_rotation()
                                ? MagneticType::type_iv
                                : MagneticType::type_iii;
  return std::make_pair(type, std::move(*rep));
}

// Transform magnetic operations by (tmat, shift) without de-duplicating:
//   (W, w) -> (tmat, shift) (W, w) (tmat, shift)^-1
// i.e. W' = tmat W tmat^-1, w' = wrap(shift - W' shift + tmat w).
[[nodiscard]] MagOps distinct_changed_magnetic_symmetry(
    Matrix3d const &tmat, Vector3d const &shift,
    std::span<MagneticSymmetryOperation const> sym_msg) {
  Matrix3d const tmat_inv = tmat.inverse();
  MagOps changed;
  changed.reserve(sym_msg.size());
  std::ranges::transform(
      sym_msg, std::back_inserter(changed), [&](auto const &op) {
        Matrix3i const rot = math::round_to_int(
            tmat * op.spatial.rotation.template cast<double>() * tmat_inv);
        Vector3d const trans = math::wrap_to_unit_cell(
            Vector3d(shift - rot.cast<double>() * shift +
                     tmat * op.spatial.translation));
        return MagneticSymmetryOperation{{rot, trans}, op.time_reversal};
      });
  return changed;
}

// Pure translations in the transformed setting: (I, w) = (tmat, shift)^-1 (I,
// w_std)(tmat, shift), i.e. w_std = tmat w.
[[nodiscard]] std::optional<std::vector<Vector3d>>
changed_pure_translations(Matrix3d const &tmat,
                          std::vector<Vector3d> const &pure_trans,
                          double symprec) {
  double const det = tmat.determinant();
  int const expected =
      math::nint(static_cast<double>(pure_trans.size()) / std::abs(det));

  std::vector<Vector3d> out;
  if (std::abs(det - 1.0) <= symprec) {
    std::ranges::transform(
        pure_trans, std::back_inserter(out),
        [&](auto const &t) { return math::wrap_to_unit_cell(Vector3d(tmat * t)); });
  } else {
    // det(tmat) need not be integer; find the least common denominator of the
    // matrix entries, then enumerate added lattice points. When no denominator
    // up to kMaxDenominator clears the entries, fail like any other unusable
    // transformation instead of silently enumerating an enormous lattice-point
    // set (the old loop fell through with denominator = kMaxDenominator + 1).
    auto const denominators = std::views::iota(1, kMaxDenominator + 1);
    auto const denom_it =
        std::ranges::find_if(denominators, [&](int denominator) {
          return std::ranges::all_of(
              std::views::cartesian_product(std::views::iota(0, 3),
                                            std::views::iota(0, 3)),
              [&](auto const st) {
                auto const [s, t] = st;
                double const scaled = tmat(s, t) * denominator;
                return std::abs(scaled - math::nint(scaled)) <= symprec;
              });
        });
    if (denom_it == denominators.end()) {
      return std::nullopt;
    }
    auto const lattice_pts = std::views::iota(0, *denom_it + 1);
    for (auto const &[n0, n1, n2, t] : std::views::cartesian_product(
             lattice_pts, lattice_pts, lattice_pts, pure_trans)) {
      Vector3d const shifted = t + Vector3d(n0, n1, n2);
      Vector3d const transformed = math::wrap_to_unit_cell(Vector3d(tmat * shifted));
      push_unique(out, transformed, [&](Vector3d const &a, Vector3d const &b) {
        return approx_equal(a, b, symprec);
      });
    }
  }

  if (static_cast<int>(out.size()) != expected) {
    return std::nullopt;
  }
  return out;
}

// The MSG operations in the changed (reference) setting, built from the coset
// representatives, the conventional centerings, and the factor group of XSG.
[[nodiscard]] std::optional<MagOps> changed_magnetic_symmetry(
    Matrix3d const &tmat, Vector3d const &shift,
    std::span<MagneticSymmetryOperation const> representatives,
    Operations const &sym_xsg, MagneticOperations const &magnetic_symmetry,
    double symprec) {
  MagOps const changed_representatives =
      distinct_changed_magnetic_symmetry(tmat, shift, representatives);

  std::vector<Vector3d> const pure_trans =
      magnetic_symmetry.pure_translations();
  auto const changed_pure =
      changed_pure_translations(tmat, pure_trans, symprec);
  if (!changed_pure) {
    return std::nullopt;
  }

  // Factor group of XSG in the conventional lattice: one operation per distinct
  // rotation, first occurrence winning.
  auto const factors = unique_by_rotation(
      sym_xsg | std::views::transform([](auto const &op) {
        return MagneticSymmetryOperation{op, false};
      }),
      [](MagneticSymmetryOperation const &op) { return op.spatial.rotation; });
  MagOps const changed_factors =
      distinct_changed_magnetic_symmetry(tmat, shift, factors);

  MagOps changed;
  changed.reserve(changed_representatives.size() * changed_pure->size() *
                  changed_factors.size());
  for (auto const &pure : *changed_pure) {
    for (auto const &rep : changed_representatives) {
      for (auto const &factor : changed_factors) {
        // (I, ti)(Pj, tj)(Pk, tk) = (Pj Pk, Pj tk + tj + ti).
        Matrix3i const rot = rep.spatial.rotation * factor.spatial.rotation;
        Vector3d const trans = math::wrap_to_unit_cell(Vector3d(
            rep.spatial.rotation.cast<double>() * factor.spatial.translation +
            rep.spatial.translation + pure));
        changed.push_back(
            {{rot, trans}, rep.time_reversal != factor.time_reversal});
      }
    }
  }
  return changed;
}

// Whether two magnetic symmetries are equal as sets (rotation, translation
// mod 1, time reversal). Rotation and time reversal select a bucket of b
// exactly; only the translations inside it are compared with tolerance.
[[nodiscard]] bool same_magnetic_symmetry(MagneticOperations const &a,
                                          MagneticOperations const &b,
                                          double symprec) {
  if (a.size() != b.size()) {
    return false;
  }
  OperationMultimap const b_index = index_by_operation_key(b.span());
  return std::ranges::all_of(a, [&](MagneticSymmetryOperation const &oa) {
    auto const [lo, hi] =
        b_index.equal_range({oa.spatial.rotation, oa.time_reversal});
    return std::ranges::any_of(lo, hi, [&](auto const &entry) {
      auto const &ob = b[static_cast<std::size_t>(entry.second)];
      // NB: strictly below symprec, unlike same_operation's <=; kept as is.
      return approx_equal(
          math::nearest_offset(
              Vector3d(oa.spatial.translation - ob.spatial.translation)),
          Vector3d::Zero(), symprec);
    });
  });
}

// Rigid rotation to the idealized standardized lattice. NB: the second factor
// is formed from `tmat^-1` (not `tmat_bravais^-1`), so the conventional-lattice
// term cancels and the result is the identity up to numerical noise; reproduced
// verbatim so the (later) std_rotation oracle matches spglib v2.7.0 exactly.
[[nodiscard]] Matrix3d rigid_rotation(Matrix3d const &lattice,
                                      Matrix3d const &tmat,
                                      SpacegroupMatch const &ref_sg) {
  Lattice const ideal_latt = refine::conventional_lattice(ref_sg);
  static_cast<void>(ideal_latt); // unused in the product (see note above)
  return lattice * tmat.inverse() * tmat * lattice.inverse();
}

struct ReferenceSpaceGroup {
  MagneticType type;
  SpacegroupMatch ref_sg;
  MagOps changed_symmetry;
  Matrix3d tmat;
  Vector3d shift;
};

// Identify FSG + XSG, the MSG type, the reference setting (FSG for types I–III,
// XSG for type IV), and the magnetic symmetry transformed into that setting.
[[nodiscard]] std::optional<ReferenceSpaceGroup>
get_reference_space_group(Matrix3d const &lattice,
                          MagneticOperations const &magnetic_symmetry,
                          double symprec) {
  auto fsg = space_group_of_magnetic_symmetry(magnetic_symmetry,
                                              SpaceGroupKind::family, symprec);
  if (!fsg) {
    return std::nullopt;
  }
  auto xsg = space_group_of_magnetic_symmetry(magnetic_symmetry,
                                              SpaceGroupKind::maximal, symprec);
  if (!xsg) {
    return std::nullopt;
  }

  auto type_rep = magnetic_space_group_type(
      magnetic_symmetry, fsg->first.size(), xsg->first.size());
  if (!type_rep) {
    return std::nullopt;
  }
  MagneticType const type = type_rep->first;
  MagOps const &representatives = type_rep->second;

  SpacegroupMatch ref_sg =
      (type == MagneticType::type_iv) ? xsg->second : fsg->second;
  Matrix3d const lattice_inv = lattice.inverse();
  ref_sg.bravais_lattice = lattice * ref_sg.bravais_lattice;
  ref_sg = refine::find_similar_bravais_lattice(ref_sg, symprec);
  ref_sg.bravais_lattice = lattice_inv * ref_sg.bravais_lattice;

  Matrix3d const tmat = ref_sg.bravais_lattice.inverse();
  Vector3d const shift = ref_sg.origin_shift;

  auto changed = changed_magnetic_symmetry(
      tmat, shift, representatives, xsg->first, magnetic_symmetry, symprec);
  if (!changed) {
    return std::nullopt;
  }
  return ReferenceSpaceGroup{type, std::move(ref_sg), std::move(*changed), tmat,
                             shift};
}

} // namespace

Result<MagneticTypeIdentification> MagneticIdentification::identify() const {
  Matrix3d const &lattice = lattice_.matrix();
  MagneticOperations const &magnetic_symmetry = operations_;
  double const symprec = tol_.symprec;
  auto reference =
      get_reference_space_group(lattice, magnetic_symmetry, symprec);
  if (!reference) {
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }
  HallNumber const hall = reference->ref_sg.hall;
  int const type = static_cast<int>(reference->type);
  auto const [first_uni, last_uni] = data::uni_candidates(hall);

  // The first UNI candidate whose tabulated operations match the changed
  // symmetry, with the correction transformation x_uni = (tmat, shift)
  // x_changed that made them match.
  struct Correction {
    UniNumber uni;
    Matrix3d tmat;
    Vector3d shift;
  };
  auto const find_correction = [&]() -> std::optional<Correction> {
    for (int u = first_uni.value(); u <= last_uni.value(); ++u) {
      UniNumber const uni = *UniNumber::of(u);
      if (data::magnetic_spacegroup_type(uni).type != type) {
        continue;
      }
      auto const &msg_uni =
          data::magnetic_operations_from_database(uni, hall);
      if (msg_uni.size() != reference->changed_symmetry.size()) {
        continue;
      }

      auto const &transformations =
          data::magnetic_std_transformations(uni, hall);
      for (auto const &transform : transformations) {
        Matrix3d const cor = transform.rotation.cast<double>();
        Vector3d const cor_shift = transform.translation;
        auto const symmetry_cor = distinct_changed_magnetic_symmetry(
            cor, cor_shift, reference->changed_symmetry);
        if (same_magnetic_symmetry(msg_uni,
                                   MagneticOperations{symmetry_cor}, symprec)) {
          return Correction{uni, cor, cor_shift};
        }
      }
    }
    return std::nullopt;
  };
  std::optional<Correction> const correction = find_correction();

  if (!correction) {
    return leaf::new_error(e_magnetic_symmetry_search_failed{});
  }

  auto const msgtype = data::magnetic_spacegroup_type(correction->uni);

  // Compose the correction onto the reference transformation:
  //   (tmat, shift) -> (tmat_cor, shift_cor) (tmat, shift).
  Matrix3d const tmat = correction->tmat * reference->tmat;
  Vector3d const shift = correction->tmat * reference->shift + correction->shift;

  SpacegroupMatch ref_sg = reference->ref_sg;
  ref_sg.bravais_lattice = lattice * ref_sg.bravais_lattice;
  Matrix3d const rigid_rot = rigid_rotation(lattice, tmat, ref_sg);

  return MagneticTypeIdentification{correction->uni,
                                    static_cast<MagneticType>(msgtype.type),
                                    hall,
                                    tmat,
                                    shift,
                                    rigid_rot};
}

Result<MagneticCell> MagneticIdentification::transform(
    MagneticCell const &mcell,
    MagneticTypeIdentification const &identification) const {
  Matrix3d const &transformation_matrix = identification.transformation_matrix;
  Vector3d const &origin_shift = identification.origin_shift;
  Matrix3d const &rigid_rotation = identification.std_rotation_matrix;
  MagneticOperations const &magnetic_symmetry = operations_;
  double const symprec = tol_.symprec;
  Cell const &cell = mcell.cell();
  bool const collinear = mcell.rank() == SiteTensor::collinear;

  // 1. Transform the cell to primitive using the magnetic pure translations.
  std::vector<Vector3d> const pure_trans =
      magnetic_symmetry.pure_translations();
  symmetry::PrimitiveFinder<GroupFamily::space> const finder(
      cell, {symprec, std::nullopt});
  BOOST_LEAF_AUTO(prim, finder.from_pure_translations(pure_trans));
  Cell const &prim_cell = prim.cell;
  // tmat_prm = tmat . cell.lattice^-1 . primitive.lattice.
  Matrix3d const tmat_prm =
      transformation_matrix * cell.lattice().matrix().inverse() * prim_cell.lattice().matrix();

  // The cell -> primitive map is many-to-one; pick the first preimage per
  // primitive atom (site tensors are invariant under pure translations).
  std::vector<std::optional<int>> remapping(
      static_cast<std::size_t>(prim_cell.size()));
  for (Index i = 0; i < cell.size(); ++i) {
    int const prim_atom = prim.mapping_table[static_cast<std::size_t>(i)];
    auto &preimage = remapping[static_cast<std::size_t>(prim_atom)];
    if (!preimage) {
      preimage = static_cast<int>(i);
    }
  }

  // 2. Pure translations of the transformed cell (the lattice points the single
  // zero translation expands into).
  auto const changed_pure =
      changed_pure_translations(tmat_prm, {Vector3d::Zero()}, symprec);
  if (!changed_pure) {
    return leaf::new_error(e_cell_standardization_failed{});
  }

  // 3. Apply the transformation and replicate across the pure translations.
  auto const np = static_cast<std::size_t>(prim_cell.size());
  std::size_t const nc = changed_pure->size();
  auto const total = static_cast<Index>(np * nc);
  Positions positions(total, 3);
  Types types(static_cast<std::size_t>(total));
  CollinearTensors scalars;
  NoncollinearTensors vectors;
  if (collinear) {
    scalars.resize(static_cast<std::size_t>(total));
  } else {
    vectors.resize(total, 3);
  }

  for (Index i = 0; i < static_cast<Index>(np); ++i) {
    Vector3d const pos_std =
        tmat_prm * prim_cell.position(i) + origin_shift; // x_std
    int const orig = *remapping[static_cast<std::size_t>(i)];
    for (std::size_t p = 0; p < nc; ++p) {
      Index const ip = i * static_cast<Index>(nc) + static_cast<Index>(p);
      auto const uip = static_cast<std::size_t>(ip);
      types[uip] = prim_cell.type(i);
      positions.row(ip) =
          math::wrap_to_unit_cell(Vector3d(pos_std + (*changed_pure)[p])).transpose();
      if (collinear) {
        scalars[uip] = mcell.scalar(orig);
      } else {
        vectors.row(ip) = (rigid_rotation * mcell.vector(orig)).transpose();
      }
    }
  }

  Lattice const lattice{rigid_rotation * cell.lattice().matrix() *
                        transformation_matrix.inverse()};
  SiteTensors tensors = collinear ? SiteTensors{std::move(scalars)}
                                  : SiteTensors{std::move(vectors)};
  return MagneticCell(Cell(lattice, positions, types), std::move(tensors));
}

} // namespace cppcrystal::magnetic
