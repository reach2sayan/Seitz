#include <spglib/spacegroup/hall_symbol.hpp>

#include <spglib/core/overlap.hpp>
#include <spglib/data/hall_generators_view.hpp>
#include <spglib/math/fractional.hpp>

#include <array>

// Port of hall_symbol.c (3D space-group path). Given the symmetry operations of
// the conventional (bravais) cell and a candidate Hall number, this determines
// whether the operations match that Hall setting and, if so, the origin shift
// that aligns them with the database operations. The Grosse-Kunstleve (1999)
// origin-shift formula shift = VSpU . dw is precomputed per setting (the VSpU
// tables); dw is the per-generator translation difference vs the database.
namespace spglib::spacegroup {

using data::Centering;

namespace {

// Centering change-of-basis matrices (hall_symbol.c) as lookup tables indexed
// by the Centering enum. M (integer) and its inverse M_inv (note M_inv^-1 == M,
// which lets transform_rotation avoid a matrix inverse). PRIMITIVE and any
// unused entry are the identity.
constexpr std::size_t kCenteringCount = 9; // Centering::error .. r_center

[[nodiscard]] Matrix3i const &centering_matrix(Centering c) {
  static auto const table = [] {
    std::array<Matrix3i, kCenteringCount> t;
    t.fill(Matrix3i::Identity());
    auto at = [&](Centering cc) -> Matrix3i & {
      return t[static_cast<std::size_t>(cc)];
    };
    at(Centering::body) << 0, 1, 1, 1, 0, 1, 1, 1, 0;
    at(Centering::face) << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    at(Centering::a_face) << 1, 0, 0, 0, 1, 1, 0, -1, 1;
    at(Centering::b_face) << 1, 0, 1, 0, 1, 0, -1, 0, 1;
    at(Centering::c_face) << 1, -1, 0, 1, 1, 0, 0, 0, 1;
    at(Centering::r_center) << 1, 0, 1, -1, 1, 1, 0, -1, 1;
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

[[nodiscard]] Matrix3d const &centering_matrix_inv(Centering c) {
  static auto const table = [] {
    std::array<Matrix3d, kCenteringCount> t;
    t.fill(Matrix3d::Identity());
    auto at = [&](Centering cc) -> Matrix3d & {
      return t[static_cast<std::size_t>(cc)];
    };
    at(Centering::body) << -0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, -0.5;
    at(Centering::face) << 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0;
    at(Centering::a_face) << 1.0, 0.0, 0.0, 0.0, 0.5, -0.5, 0.0, 0.5, 0.5;
    at(Centering::b_face) << 0.5, 0.0, -0.5, 0.0, 1.0, 0.0, 0.5, 0.0, 0.5;
    at(Centering::c_face) << 0.5, 0.5, 0.0, -0.5, 0.5, 0.0, 0.0, 0.0, 1.0;
    at(Centering::r_center) << 2. / 3, -1. / 3, -1. / 3, 1. / 3, 1. / 3,
        -2. / 3, 1. / 3, 1. / 3, 1. / 3;
    return t;
  }();
  return table[static_cast<std::size_t>(c)];
}

// trans expressed in the primitive setting: M . trans (transform_translation).
[[nodiscard]] Vector3d transform_translation(Centering c,
                                             Vector3d const &trans) {
  return centering_matrix(c).cast<double>() * trans;
}

// rot expressed in the primitive setting: M . rot . M_inv (transform_rotation,
// the similarity transform b^-1 . rot . b with b = M_inv, b^-1 = M).
[[nodiscard]] Matrix3d transform_rotation(Centering c, Matrix3i const &rot) {
  return centering_matrix(c).cast<double>() * rot.cast<double>() *
         centering_matrix_inv(c);
}

[[nodiscard]] std::array<Matrix3i, 3>
unpack_generators(data::GeneratorSet const &g) {
  std::array<Matrix3i, 3> rot;
  for (std::size_t i = 0; i < 3; ++i)
    rot[i] = data::generator_matrix(g, i);
  return rot;
}

// Find, for each non-zero generator rotation, a matching operation in
// `symmetry` and record its translation (get_translations). Returns false if a
// generator rotation is absent.
[[nodiscard]] bool get_translations(std::array<Vector3d, 3> &trans,
                                    SymmetryOperations const &symmetry,
                                    std::array<Matrix3i, 3> const &rot) {
  trans = {Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero()};
  for (int i = 0; i < 3; ++i) {
    if (rot[static_cast<std::size_t>(i)].isZero())
      continue;
    bool found = false;
    for (auto const &op : symmetry)
      if (op.rotation == rot[static_cast<std::size_t>(i)]) {
        trans[static_cast<std::size_t>(i)] = op.translation;
        found = true;
        break;
      }
    if (!found)
      return false;
  }
  return true;
}

// dw for one generator: translation difference (primitive setting) vs the
// matching database operation (set_dw).
[[nodiscard]] bool set_dw(Vector3d &dw, SymmetryOperations const &db_ops,
                          Matrix3i const &rot, Vector3d const &trans,
                          Centering c) {
  Vector3d const trans_prim = transform_translation(c, trans);
  for (auto const &db : db_ops)
    if (db.rotation == rot) {
      dw = trans_prim - transform_translation(c, db.translation);
      return true;
    }
  return false;
}

// shift = VSpU . dw (get_origin_shift), with dw assembled per generator.
[[nodiscard]] bool get_origin_shift(Vector3d &shift, int hall_number,
                                    std::array<Matrix3i, 3> const &rot,
                                    std::array<Vector3d, 3> const &trans,
                                    Centering c, data::VSpUSet const &vspu) {
  auto const &db_ops = data::database().operations.at(hall_number);
  std::array<double, 9> dw{};
  for (int i = 0; i < 3; ++i) {
    auto const ui = static_cast<std::size_t>(i);
    if (rot[ui].determinant() == 0)
      continue; // zero generator -> dw block stays 0
    Vector3d tmp;
    if (!set_dw(tmp, db_ops, rot[ui], trans[ui], c))
      return false;
    // hall_number is always > 0 here (3D; layer groups, hall_number < 0, would
    // leave the c component un-folded — out of scope).
    dw[ui * 3] = math::mod1(tmp[0]);
    dw[ui * 3 + 1] = math::mod1(tmp[1]);
    dw[ui * 3 + 2] = math::mod1(tmp[2]);
  }

  shift = data::vspu_matrix(vspu) * Eigen::Map<data::DwVector const>(dw.data());
  shift = math::mod1(shift);
  return true;
}

// Verify every operation in `symmetry` matches a database operation once the
// origin shift is applied (is_match_database).
[[nodiscard]] bool is_match_database(int hall_number, Vector3d const &shift,
                                     Matrix3d const &primitive_lattice,
                                     Centering c,
                                     SymmetryOperations const &symmetry,
                                     double symprec) {
  auto const &db_ops = data::database().operations.at(hall_number);
  std::vector<char> found(db_ops.size(), 0);
  for (auto const &op : symmetry) {
    bool is_found = false;
    for (std::size_t j = 0; j < db_ops.size(); ++j) {
      if (op.rotation != db_ops[j].rotation)
        continue;
      Vector3d const diff = transform_translation(c, op.translation) -
                            transform_translation(c, db_ops[j].translation) +
                            shift;
      Vector3d const shift_rot =
          transform_rotation(c, db_ops[j].rotation) * shift;
      if (!found[j] &&
          is_overlap(diff, shift_rot, primitive_lattice, symprec)) {
        found[j] = 1;
        is_found = true;
        break;
      }
    }
    if (!is_found)
      return false;
  }
  return true;
}

// Core single-candidate test (is_hall_symbol).
[[nodiscard]] bool is_hall_symbol(Vector3d &shift, int hall_number,
                                  Matrix3d const &primitive_lattice,
                                  SymmetryOperations const &symmetry,
                                  Centering c, data::GeneratorSet const &gens,
                                  data::VSpUSet const &vspu, double symprec) {
  auto const &db_ops = data::database().operations.at(hall_number);
  if (db_ops.size() != symmetry.size())
    return false;
  auto const rot = unpack_generators(gens);
  std::array<Vector3d, 3> trans;
  if (!get_translations(trans, symmetry, rot))
    return false;
  if (!get_origin_shift(shift, hall_number, rot, trans, c, vspu))
    return false;
  return is_match_database(hall_number, shift, primitive_lattice, c, symmetry,
                           symprec);
}

// ---- per-crystal-system candidate loops ----

template <std::size_t N>
[[nodiscard]] bool
try_entries(Vector3d &shift, int hall_number, Matrix3d const &prim,
            SymmetryOperations const &sym, Centering c,
            std::array<data::GeneratorSet, N> const &gens,
            std::array<data::VSpUSet, N> const &vspu, double symprec) {
  for (std::size_t i = 0; i < N; ++i)
    if (is_hall_symbol(shift, hall_number, prim, sym, c, gens[i], vspu[i],
                       symprec))
      return true;
  return false;
}

[[nodiscard]] bool is_rhombohedral_hall(int h) {
  switch (h) {
  case 433:
  case 434:
  case 436:
  case 437:
  case 444:
  case 445:
  case 450:
  case 451:
  case 452:
  case 453:
  case 458:
  case 459:
  case 460:
  case 461:
    return true;
  default:
    return false;
  }
}

[[nodiscard]] bool is_rhombo_h_setting(int h) {
  switch (h) {
  case 433:
  case 436:
  case 444:
  case 450:
  case 452:
  case 458:
  case 460:
    return true;
  default:
    return false;
  }
}

[[nodiscard]] bool dispatch(Vector3d &shift, int hall_number,
                            Matrix3d const &prim, SymmetryOperations const &sym,
                            Centering c, double symprec) {
  using namespace data;
  if (489 <= hall_number && hall_number <= 530) { // cubic
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
  }
  if (462 <= hall_number && hall_number <= 488) // hexagonal
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       hexa_generators, hexa_VSpU, symprec);
  if (430 <= hall_number && hall_number <= 461) { // trigonal / rhombohedral
    if (is_rhombohedral_hall(hall_number)) {
      if (is_rhombo_h_setting(hall_number))
        return try_entries(shift, hall_number, prim, sym, Centering::r_center,
                           rhombo_h_generators, rhombo_h_VSpU, symprec);
      return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                         rhombo_p_generators, rhombo_p_VSpU, symprec);
    }
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       trigo_generators, trigo_VSpU, symprec);
  }
  if (349 <= hall_number && hall_number <= 429) { // tetragonal
    if (c == Centering::primitive)
      return try_entries(shift, hall_number, prim, sym, c, tetra_generators,
                         tetra_VSpU, symprec);
    if (c == Centering::body)
      return try_entries(shift, hall_number, prim, sym, c, tetra_generators,
                         tetra_I_VSpU, symprec);
    return false;
  }
  if (108 <= hall_number && hall_number <= 348) { // orthorhombic
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
  }
  if (3 <= hall_number && hall_number <= 107) { // monoclinic
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
  }
  if (1 <= hall_number && hall_number <= 2) // triclinic
    return try_entries(shift, hall_number, prim, sym, Centering::primitive,
                       tricli_generators, tricli_VSpU, symprec);
  return false;
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

} // namespace spglib::spacegroup
