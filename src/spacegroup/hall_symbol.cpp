#include <cppcrystal/spacegroup/hall_symbol.hpp>

#include <cppcrystal/core/centering.hpp>
#include <cppcrystal/core/matrix_order.hpp>
#include <cppcrystal/core/overlap.hpp>
#include <cppcrystal/core/periodicity.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/data/hall_classification.hpp>
#include <cppcrystal/data/hall_generators_view.hpp>

#include <boost/container/static_vector.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>

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

// Everything fixed for one (Hall number, operation set) match: both operation
// lists with their rotation indices, so every "the operation with rotation R"
// question below is a logarithmic lookup.
struct Setting {
  int hall_number;
  Matrix3d const &primitive_lattice;
  SymmetryOperations const &symmetry;
  RotationMultimap<int> const &symmetry_by_rotation;
  SymmetryOperations const &db_ops;
  RotationMultimap<int> const &db_by_rotation;
  double symprec;

  // For a layer setting (hall_number < 0) the conventional c axis (index 2)
  // is aperiodic: the overlap test must not fold it, and a shift along it
  // stays raw — that is what re-centers a layer onto the database convention.
  [[nodiscard]] std::optional<int> aperiodic_axis() const noexcept {
    return hall_number < 0 ? std::optional<int>(2) : std::nullopt;
  }
  [[nodiscard]] Vector3d wrap(Vector3d const &v) const noexcept {
    return wrap_periodic(v, aperiodic_axis());
  }
};

using Generators = std::array<Matrix3i, 3>;
using Translations = std::array<Vector3d, 3>;

[[nodiscard]] Generators unpack_generators(data::GeneratorSet const &g) {
  return {data::generator_matrix(g, 0), data::generator_matrix(g, 1),
          data::generator_matrix(g, 2)};
}

// The translation of the first operation carrying each non-zero generator
// rotation; nullopt if a generator rotation is absent from the operations.
[[nodiscard]] std::optional<Translations>
generator_translations(Setting const &s, Generators const &rot) {
  Translations trans{Vector3d::Zero(), Vector3d::Zero(), Vector3d::Zero()};
  for (auto [t, r] : std::views::zip(trans, rot)) {
    if (r.isZero()) {
      continue;
    }
    auto const it = s.symmetry_by_rotation.find(r);
    if (it == s.symmetry_by_rotation.end()) {
      return std::nullopt;
    }
    t = s.symmetry[static_cast<std::size_t>(it->second)].translation;
  }
  return trans;
}

// dw for one generator: translation difference (primitive setting) vs the
// first database operation with that rotation.
[[nodiscard]] std::optional<Vector3d> dw_of(Setting const &s, Centering c,
                                            Matrix3i const &rot,
                                            Vector3d const &trans) {
  auto const it = s.db_by_rotation.find(rot);
  if (it == s.db_by_rotation.end()) {
    return std::nullopt;
  }
  auto const &db = s.db_ops[static_cast<std::size_t>(it->second)];
  return Vector3d(transform_translation(c, trans) -
                  transform_translation(c, db.translation));
}

// shift = VSpU . dw, with dw assembled per generator (a zero generator
// contributes a zero block).
[[nodiscard]] std::optional<Vector3d>
origin_shift(Setting const &s, Centering c, Generators const &rot,
             Translations const &trans, data::VSpUSet const &vspu) {
  data::DwVector dw = data::DwVector::Zero();
  for (auto const [i, r, t] : std::views::zip(std::views::iota(0, 3), rot, trans)) {
    if (r.determinant() == 0) {
      continue;
    }
    auto const d = dw_of(s, c, r, t);
    if (!d) {
      return std::nullopt;
    }
    dw.segment<3>(static_cast<Index>(3 * i)) = s.wrap(*d);
  }
  return s.wrap(data::vspu_matrix(vspu) * dw);
}

// Every operation must reproduce a distinct database operation once the
// origin shift is applied: same rotation, translations agreeing within
// symprec. Operations are matched in order, each to the first still-unmatched
// database operation carrying its rotation.
[[nodiscard]] bool matches_database(Setting const &s, Centering c,
                                    Vector3d const &shift) {
  boost::container::static_vector<bool, 192> matched(s.db_ops.size(), false);
  for (auto const &op : s.symmetry) {
    Vector3d const lhs = transform_translation(c, op.translation) + shift;
    Vector3d const shift_rot = transform_rotation(c, op.rotation) * shift;
    auto const [lo, hi] = s.db_by_rotation.equal_range(op.rotation);
    auto const it = std::ranges::find_if(lo, hi, [&](auto const &entry) {
      auto const idx = static_cast<std::size_t>(entry.second);
      return !matched[idx] &&
             is_overlap(lhs - transform_translation(c, s.db_ops[idx].translation),
                        shift_rot, s.primitive_lattice, s.symprec,
                        s.aperiodic_axis());
    });
    if (it == hi) {
      return false;
    }
    matched[static_cast<std::size_t>(it->second)] = true;
  }
  return true;
}

// One (generators, VSpU) candidate: the origin shift if the operations match
// the setting through it.
[[nodiscard]] std::optional<Vector3d>
hall_symbol_shift(Setting const &s, Centering c, data::GeneratorSet const &gens,
                  data::VSpUSet const &vspu) {
  Generators const rot = unpack_generators(gens);
  return generator_translations(s, rot)
      .and_then([&](Translations const &trans) {
        return origin_shift(s, c, rot, trans, vspu);
      })
      .and_then([&](Vector3d const &shift) {
        return matches_database(s, c, shift) ? std::optional<Vector3d>(shift)
                                             : std::nullopt;
      });
}

// The first candidate of a family that matches.
[[nodiscard]] std::optional<Vector3d>
first_shift(Setting const &s, Centering c,
            std::span<data::GeneratorSet const> gens,
            std::span<data::VSpUSet const> vspu) {
  for (auto const &[gen, vsp] : std::views::zip(gens, vspu)) {
    if (auto shift = hall_symbol_shift(s, c, gen, vsp)) {
      return shift;
    }
  }
  return std::nullopt;
}

// The Grosse-Kunstleve candidate family for one (holohedry, centering) cell:
// the generator and VSpU tables plus the centering the matching runs in.
struct Family {
  std::span<data::GeneratorSet const> gens;
  std::span<data::VSpUSet const> vspu;
  Centering pass;
};

constexpr std::size_t kNumHolohedries = 8; // Holohedry::none .. cubic
constexpr std::size_t kNumCenterings = 9;  // Centering::error .. r_center

// (holohedry, input centering) -> family, or nullopt where the tables have no
// candidate. Rows keyed on the holohedry alone (systems whose matching always
// runs in the primitive setting) fill every centering column.
constexpr auto kFamilies = [] {
  using namespace data;
  using enum Centering;
  std::array<std::array<std::optional<Family>, kNumCenterings>, kNumHolohedries>
      t{};
  auto const row = [&](Holohedry h, std::optional<Centering> match,
                       std::span<GeneratorSet const> gens,
                       std::span<VSpUSet const> vspu,
                       std::optional<Centering> pass) {
    for (std::size_t c = 0; c < kNumCenterings; ++c) {
      auto const centering = static_cast<Centering>(c);
      if (!match || *match == centering) {
        t[static_cast<std::size_t>(h)][c] =
            Family{gens, vspu, pass.value_or(centering)};
      }
    }
  };
  constexpr std::optional<Centering> any = std::nullopt;
  constexpr std::optional<Centering> keep = std::nullopt;
  row(Holohedry::cubic, primitive, cubic_generators, cubic_VSpU, keep);
  row(Holohedry::cubic, body, cubic_generators, cubic_I_VSpU, keep);
  row(Holohedry::cubic, face, cubic_generators, cubic_F_VSpU, keep);
  row(Holohedry::hexagonal, any, hexa_generators, hexa_VSpU, primitive);
  row(Holohedry::trigonal, any, trigo_generators, trigo_VSpU, primitive);
  row(Holohedry::tetragonal, primitive, tetra_generators, tetra_VSpU, keep);
  row(Holohedry::tetragonal, body, tetra_generators, tetra_I_VSpU, keep);
  row(Holohedry::orthorhombic, primitive, ortho_generators, ortho_VSpU, keep);
  row(Holohedry::orthorhombic, body, ortho_generators, ortho_I_VSpU, keep);
  row(Holohedry::orthorhombic, face, ortho_generators, ortho_F_VSpU, keep);
  row(Holohedry::orthorhombic, a_face, ortho_generators, ortho_A_VSpU, keep);
  row(Holohedry::orthorhombic, b_face, ortho_generators, ortho_B_VSpU, keep);
  row(Holohedry::orthorhombic, c_face, ortho_generators, ortho_C_VSpU, keep);
  row(Holohedry::monoclinic, primitive, monocli_generators, monocli_VSpU, keep);
  row(Holohedry::monoclinic, a_face, monocli_generators, monocli_A_VSpU, keep);
  row(Holohedry::monoclinic, b_face, monocli_generators, monocli_B_VSpU, keep);
  row(Holohedry::monoclinic, c_face, monocli_generators, monocli_C_VSpU, keep);
  row(Holohedry::monoclinic, body, monocli_generators, monocli_I_VSpU, keep);
  row(Holohedry::triclinic, any, tricli_generators, tricli_VSpU, primitive);
  return t;
}();

[[nodiscard]] std::optional<Vector3d> dispatch(Setting const &s, Centering c) {
  using namespace data;
  // Crystal system and the rhombohedral subsets are pure functions of the Hall
  // number, precomputed once in data::kHallClass. Layer groups (negative hall
  // numbers) carry no 3D hall-range bucket, so their crystal system is derived
  // from the point-group number; they are never rhombohedral and reuse the
  // same per-system generator families.
  HallClass const hc =
      s.hall_number < 0
          ? HallClass{holohedry_from_pointgroup(
                          spacegroup_type(s.hall_number).pointgroup_number),
                      false, false}
          : hall_class(s.hall_number);

  // The one choice not keyed on (system, centering): rhombohedral settings
  // split on the hexagonal-vs-primitive axes choice of the Hall symbol.
  if (hc.system == Holohedry::trigonal && hc.rhombohedral) {
    return hc.rhombo_hex_setting
               ? first_shift(s, Centering::r_center, rhombo_h_generators,
                             rhombo_h_VSpU)
               : first_shift(s, Centering::primitive, rhombo_p_generators,
                             rhombo_p_VSpU);
  }

  auto const &family = kFamilies[static_cast<std::size_t>(hc.system)]
                                [static_cast<std::size_t>(c)];
  if (!family) {
    return std::nullopt;
  }
  return first_shift(s, family->pass, family->gens, family->vspu);
}

} // namespace

std::optional<Vector3d> match_hall_symbol(Matrix3d const &bravais_lattice,
                                          int hall_number, Centering centering,
                                          SymmetryOperations const &symmetry,
                                          double symprec) {
  SymmetryOperations const &db_ops = data::operations_from_database(hall_number);
  if (db_ops.size() != symmetry.size()) {
    return std::nullopt;
  }
  Matrix3d const primitive_lattice =
      bravais_lattice * centering_matrix_inv(centering);
  RotationMultimap<int> const symmetry_by_rotation =
      index_by_rotation(symmetry, &SymmetryOperation::rotation);
  Setting const setting{hall_number,
                        primitive_lattice,
                        symmetry,
                        symmetry_by_rotation,
                        db_ops,
                        data::operations_by_rotation(hall_number),
                        symprec};
  // Bring the origin shift back to the bravais setting.
  return dispatch(setting, centering).transform([&](Vector3d const &shift) {
    return Vector3d(centering_matrix_inv(centering) * shift);
  });
}

} // namespace cppcrystal::spacegroup
