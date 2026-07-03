#include <cppcrystal/spacegroup/hall_symbol.hpp>

#include <cppcrystal/core/centering.hpp>
#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/data/hall_classification.hpp>
#include <cppcrystal/data/hall_generators_view.hpp>
#include <cppcrystal/math/fractional.hpp>

#include <boost/container/static_vector.hpp>

#include <array>
#include <cstdint>
#include <ranges>

// Given the symmetry operations of the conventional (bravais) cell and a
// candidate Hall number, determine whether the operations match that Hall
// setting and, if so, the origin shift that aligns them with the database
// operations. The Grosse-Kunstleve (1999) origin-shift formula shift = VSpU . dw
// is precomputed per setting (the VSpU tables); dw is the per-generator
// translation difference vs the database.
namespace cppcrystal::spacegroup {

using data::Centering;

namespace {

// Centering change-of-basis matrices live in core/centering.hpp (shared with the
// space-group search and cell standardization).

// trans expressed in the primitive setting: M . trans.
[[nodiscard]] Vector3d transform_translation(Centering c,
                                             Vector3d const &trans) {
  return centering_matrix(c).cast<double>() * trans;
}

// rot expressed in the primitive setting: M . rot . M_inv (the similarity
// transform b^-1 . rot . b with b = M_inv, b^-1 = M).
[[nodiscard]] Matrix3d transform_rotation(Centering c, Matrix3i const &rot) {
  return centering_matrix(c).cast<double>() * rot.cast<double>() *
         centering_matrix_inv(c);
}

[[nodiscard]] std::array<Matrix3i, 3>
unpack_generators(data::GeneratorSet const &g) {
  std::array<Matrix3i, 3> rot;
  for (std::size_t i = 0; i < 3; ++i) {
    rot[i] = data::generator_matrix(g, i);
  }
  return rot;
}

// Find, for each non-zero generator rotation, a matching operation in
// `symmetry` and record its translation. Returns false if a generator rotation
// is absent.
[[nodiscard]] bool get_translations(std::array<Vector3d, 3> &trans,
                                    SymmetryOperations const &symmetry,
                                    std::array<Matrix3i, 3> const &rot) {
  trans = {Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero()};

  for (auto [t, r] : std::views::zip(trans, rot)) {
    if (r.isZero()) {
      continue;
    }

    const auto it = std::ranges::find_if(
        symmetry, [&](const auto &op) { return op.rotation == r; });
    if (it == symmetry.end()) {
      return false;
    }
    t = it->translation;
  }

  return true;
}

// dw for one generator: translation difference (primitive setting) vs the
// matching database operation.
[[nodiscard]] bool set_dw(Vector3d &dw, SymmetryOperations const &db_ops,
                          Matrix3i const &rot, Vector3d const &trans,
                          Centering c) {
  const Vector3d trans_prim = transform_translation(c, trans);
  const auto it = std::ranges::find_if(
      db_ops, [&](const auto &db) { return db.rotation == rot; });
  if (it == db_ops.end()) {
    return false;
  }

  dw = trans_prim - transform_translation(c, it->translation);
  return true;
}

// Fold a fractional 3-vector into the cell. For a layer setting (hall_number <
// 0) the conventional c axis (index 2) is the aperiodic axis and is not
// periodic, so its component is left raw rather than folded into [0, 1) — this
// is what lets a layer offset along its aperiodic axis be re-centered onto the
// database convention (the third component is `hall_number > 0 ? folded : raw`).
[[nodiscard]] Vector3d wrap_shift(Vector3d const &v, int hall_number) {
  Vector3d out = math::wrap_to_unit_cell(v);
  if (hall_number < 0) {
    out[2] = v[2];
  }
  return out;
}

// shift = VSpU . dw, with dw assembled per generator.
[[nodiscard]] bool get_origin_shift(Vector3d &shift, int hall_number,
                                    std::array<Matrix3i, 3> const &rot,
                                    std::array<Vector3d, 3> const &trans,
                                    Centering c, data::VSpUSet const &vspu) {
  auto const &db_ops = data::operations_from_database(hall_number);
  data::DwVector dw = data::DwVector::Zero();
  for (std::size_t i = 0; i < 3; ++i) {
    if (rot[i].determinant() == 0) {
      continue; // zero generator -> dw block stays 0
    }
    Vector3d tmp;
    if (!set_dw(tmp, db_ops, rot[i], trans[i], c)) {
      return false;
    }

    dw.segment<3>(static_cast<Eigen::Index>(i * 3)) =
        wrap_shift(tmp, hall_number);
  }

  shift = wrap_shift(data::vspu_matrix(vspu) * dw, hall_number);
  return true;
}

// Verify every operation in `symmetry` matches a database operation once the
// origin shift is applied.
[[nodiscard]] bool is_match_database(int hall_number, Vector3d const &shift,
                                     Matrix3d const &primitive_lattice,
                                     Centering c,
                                     SymmetryOperations const &symmetry,
                                     double symprec) {
  auto const &db_ops = data::operations_from_database(hall_number);
  // For a layer setting the aperiodic axis (conventional c = 2) is not periodic,
  // so the overlap test must not fold it.
  std::optional<int> const aperiodic_axis =
      hall_number < 0 ? std::optional<int>(2) : std::nullopt;
  boost::container::static_vector<std::uint8_t, 192> found(db_ops.size(), 0);
  for (const auto &op : symmetry) {
    const auto it = std::ranges::find_if(
        db_ops | std::views::enumerate, [&](const auto &entry) {
          const auto [idx, db_op] = entry;

          if (found[static_cast<std::size_t>(idx)] ||
              op.rotation != db_op.rotation) {
            return false;
          }

          const Vector3d diff = transform_translation(c, op.translation) -
                                transform_translation(c, db_op.translation) +
                                shift;

          const Vector3d shift_rot =
              transform_rotation(c, db_op.rotation) * shift;

          return is_overlap(diff, shift_rot, primitive_lattice, symprec,
                            aperiodic_axis);
        });

    if (it == std::ranges::end(db_ops | std::views::enumerate)) {
      return false;
    }

    found[static_cast<std::size_t>(std::get<0>(*it))] = 1;
  }

  return true;
}

// Core single-candidate test.
[[nodiscard]] bool is_hall_symbol(Vector3d &shift, int hall_number,
                                  Matrix3d const &primitive_lattice,
                                  SymmetryOperations const &symmetry,
                                  Centering c, data::GeneratorSet const &gens,
                                  data::VSpUSet const &vspu, double symprec) {
  auto const &db_ops = data::operations_from_database(hall_number);
  if (db_ops.size() != symmetry.size()) {
    return false;
  }
  auto const rot = unpack_generators(gens);
  std::array<Vector3d, 3> trans;

  bool found = get_translations(trans, symmetry, rot);
  if (found) {
    found = get_origin_shift(shift, hall_number, rot, trans, c, vspu);
  }
  if (found) {
    found = is_match_database(hall_number, shift, primitive_lattice, c,
                              symmetry, symprec);
  }
  return found;
}

// ---- per-crystal-system candidate loops ----
template <std::size_t N>
[[nodiscard]] bool
try_entries(Vector3d &shift, int hall_number, Matrix3d const &prim,
            SymmetryOperations const &sym, Centering c,
            std::array<data::GeneratorSet, N> const &gens,
            std::array<data::VSpUSet, N> const &vspu, double symprec) {
  return std::ranges::any_of(std::views::zip(gens, vspu), [&](const auto &x) {
    const auto &[gen, vsp] = x;
    return is_hall_symbol(shift, hall_number, prim, sym, c, gen, vsp, symprec);
  });
}

[[nodiscard]] bool dispatch(Vector3d &shift, int hall_number,
                            Matrix3d const &prim, SymmetryOperations const &sym,
                            Centering c, double symprec) {
  using namespace data;
  // Crystal system and the rhombohedral subsets are pure functions of the Hall
  // number, precomputed once in data::kHallClass. The VSpU/generator choice
  // below still depends on the runtime Centering c, so it stays here. Layer
  // groups (negative hall numbers) carry no 3D hall-range bucket, so their
  // crystal system is derived from the point-group number; they are never
  // rhombohedral and reuse the same per-system generator families.
  HallClass const hc =
      hall_number < 0
          ? HallClass{holohedry_from_pointgroup(
                          spacegroup_type(hall_number).pointgroup_number),
                      false, false}
          : hall_class(hall_number);
  switch (hc.system) {
  case Holohedry::cubic:
    if (c == Centering::primitive)
      return try_entries(shift, hall_number, prim, sym, c, cubic_generators,
                         cubic_VSpU, symprec);
    if (c == Centering::body)
      return try_entries(shift, hall_number, prim, sym, c, cubic_generators,
                         cubic_I_VSpU, symprec);
    if (c == Centering::face)
      return try_entries(shift, hall_number, prim, sym, c, cubic_generators,
                         cubic_F_VSpU, symprec);
    return false;
  case Holohedry::hexagonal:
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       hexa_generators, hexa_VSpU, symprec);
  case Holohedry::trigonal:
    if (hc.rhombohedral) {
      if (hc.rhombo_hex_setting)
        return try_entries(shift, hall_number, prim, sym, Centering::r_center,
                           rhombo_h_generators, rhombo_h_VSpU, symprec);
      return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                         rhombo_p_generators, rhombo_p_VSpU, symprec);
    }
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       trigo_generators, trigo_VSpU, symprec);
  case Holohedry::tetragonal:
    if (c == Centering::primitive)
      return try_entries(shift, hall_number, prim, sym, c, tetra_generators,
                         tetra_VSpU, symprec);
    if (c == Centering::body)
      return try_entries(shift, hall_number, prim, sym, c, tetra_generators,
                         tetra_I_VSpU, symprec);
    return false;
  case Holohedry::orthorhombic:
    switch (c) {
    case Centering::primitive:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_VSpU, symprec);
    case Centering::body:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_I_VSpU, symprec);
    case Centering::face:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_F_VSpU, symprec);
    case Centering::a_face:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_A_VSpU, symprec);
    case Centering::b_face:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_B_VSpU, symprec);
    case Centering::c_face:
      return try_entries(shift, hall_number, prim, sym, c, ortho_generators,
                         ortho_C_VSpU, symprec);
    default:
      return false;
    }
  case Holohedry::monoclinic:
    switch (c) {
    case Centering::primitive:
      return try_entries(shift, hall_number, prim, sym, c, monocli_generators,
                         monocli_VSpU, symprec);
    case Centering::a_face:
      return try_entries(shift, hall_number, prim, sym, c, monocli_generators,
                         monocli_A_VSpU, symprec);
    case Centering::b_face:
      return try_entries(shift, hall_number, prim, sym, c, monocli_generators,
                         monocli_B_VSpU, symprec);
    case Centering::c_face:
      return try_entries(shift, hall_number, prim, sym, c, monocli_generators,
                         monocli_C_VSpU, symprec);
    case Centering::body:
      return try_entries(shift, hall_number, prim, sym, c, monocli_generators,
                         monocli_I_VSpU, symprec);
    default:
      return false;
    }
  case Holohedry::triclinic:
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       tricli_generators, tricli_VSpU, symprec);
  case Holohedry::none:
  default:
    return false;
  }
}

} // namespace

std::optional<Vector3d> match_hall_symbol(Matrix3d const &bravais_lattice,
                                          int hall_number, Centering centering,
                                          SymmetryOperations const &symmetry,
                                          double symprec) {
  Matrix3d const primitive_lattice =
      bravais_lattice * centering_matrix_inv(centering);
  Vector3d shift;
  if (!dispatch(shift, hall_number, primitive_lattice, symmetry, centering,
                symprec))
    return std::nullopt;
  // Bring the origin shift back to the bravais setting.
  return Vector3d(centering_matrix_inv(centering) * shift);
}

} // namespace cppcrystal::spacegroup
