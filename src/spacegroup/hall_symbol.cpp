#include "spacegroup/spacegroup.hpp"

#include "core/centering.hpp"
#include "core/matrix_order.hpp"
#include "core/position_index.hpp"
#include "data/hall_classification.hpp"
#include "data/hall_generators_view.hpp"
#include "data/operation_index.hpp"
#include <seitz/core/operation_set.hpp>
#include <seitz/core/periodicity.hpp>
#include <seitz/core/point_group.hpp>
#include <seitz/core/tolerance.hpp>

#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>

#include <array>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

// Given the symmetry operations of the conventional (bravais) cell and a
// candidate Hall number, determine whether the operations match that Hall
// setting and, if so, the origin shift that aligns them with the database
// operations. The Grosse-Kunstleve (1999) origin-shift formula shift = VSpU .
// dw is precomputed per setting (the VSpU tables); dw is the per-generator
// translation difference vs the database, over every representative of the
// generator rotation the operations offer.
namespace seitz::spacegroup {

using data::Centering;

namespace {

// Centering change-of-basis matrices live in core/centering.hpp (shared with
// the space-group search and cell standardization).

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

// Everything fixed for one (Hall setting, operation set) match: both operation
// lists with their rotation indices, so every "the operation with rotation R"
// question below is a logarithmic lookup. The family is a template parameter,
// so the layer path's aperiodic c axis is an `if constexpr` branch.
template <GroupFamily F> struct MatchContext {
  HallNumber hall;
  Matrix3d const &primitive_lattice;
  Operations const &symmetry;
  RotationMultimap<int> const &symmetry_by_rotation;
  Operations const &db_ops;
  RotationMultimap<int> const &db_by_rotation;
  double symprec;

  // For a layer setting the conventional c axis (index 2) is aperiodic: the
  // overlap test must not fold it, and a shift along it stays raw — that is
  // what re-centers a layer onto the database convention.
  [[nodiscard]] static constexpr CellPeriodicity periodicity() noexcept {
    if constexpr (F == GroupFamily::layer) {
      return aperiodic_along(2);
    } else {
      return all_periodic();
    }
  }
  [[nodiscard]] Vector3d wrap(Vector3d const &v) const noexcept {
    return seitz::wrap(v, periodicity());
  }
};

using Generators = std::array<Matrix3i, 3>;

// The dw values one generator rotation admits, in operation order, and the
// three of them a generator set needs.
using DwCandidates = boost::container::small_vector<Vector3d, 2>;
using DwChoices = std::array<DwCandidates, 3>;

[[nodiscard]] Generators unpack_generators(data::GeneratorSet const &g) {
  return {data::generator_matrix(g, 0), data::generator_matrix(g, 1),
          data::generator_matrix(g, 2)};
}

// Every distinct dw a generator rotation admits: the folded translation
// difference (primitive setting) between an operation with that rotation and
// the first database operation carrying it. Two operations of one rotation
// differ by a lattice translation, integral in the primitive setting, so a
// declared centering folds them onto one value. The B-centered monoclinic and
// orthorhombic settings (13, 15, 34, ... 333) are declared PRIMITIVE while
// their operations still carry the B translation, so each representative gives
// its own dw and only trying all of them makes the match independent of input
// order. A zero generator contributes the single zero block.
template <GroupFamily F>
[[nodiscard]] std::optional<DwCandidates>
dw_candidates(MatchContext<F> const &s, Centering c, Matrix3i const &rot) {
  if (rot.determinant() == 0) {
    return DwCandidates{Vector3d::Zero()};
  }
  auto const db = s.db_by_rotation.find(rot);
  if (db == s.db_by_rotation.end()) {
    return std::nullopt;
  }
  Vector3d const reference = transform_translation(
      c, s.db_ops[static_cast<std::size_t>(db->second)].translation);
  auto const [lo, hi] = s.symmetry_by_rotation.equal_range(rot);
  DwCandidates out;
  for (int const index : std::ranges::subrange(lo, hi) | std::views::values) {
    push_unique(
        out,
        s.wrap(transform_translation(
                   c, s.symmetry[static_cast<std::size_t>(index)].translation) -
               reference),
        [](Vector3d const &a, Vector3d const &b) {
          return approx_equal(a, b, kZeroPrec);
        });
  }
  if (out.empty()) {
    return std::nullopt; // the generator rotation is absent from `symmetry`
  }
  return out;
}

// The per-generator dw candidates; nullopt if any generator rotation is
// missing from either operation list.
template <GroupFamily F>
[[nodiscard]] std::optional<DwChoices>
dw_choices(MatchContext<F> const &s, Centering c, Generators const &rot) {
  DwChoices choices;
  for (auto [out, r] : std::views::zip(choices, rot)) {
    auto candidates = dw_candidates(s, c, r);
    if (!candidates) {
      return std::nullopt;
    }
    out = std::move(*candidates);
  }
  return choices;
}

// shift = VSpU . dw for one choice of generator representatives.
template <GroupFamily F>
[[nodiscard]] Vector3d
origin_shift(MatchContext<F> const &s, data::VSpUSet const &vspu,
             Vector3d const &dw0, Vector3d const &dw1, Vector3d const &dw2) {
  data::DwVector dw;
  dw << dw0, dw1, dw2;
  return s.wrap(data::vspu_matrix(vspu) * dw);
}

// The two operation lists re-expressed in the primitive setting: everything
// matches_database needs that does NOT depend on the origin shift. It runs once
// per element of the dw product below, and rebuilding M.t and M.R.M^-1 for up
// to 192 operations each time was the bulk of its cost. Built once per
// (operation set, centering), i.e. once per match_hall.
struct CenteredOperations {
  std::vector<Vector3d> sym_translation; // M . t, per input operation
  std::vector<Matrix3d> sym_rotation;    // M . R . M^-1, per input operation
  std::vector<Vector3d> db_translation;  // M . t, per database operation
};

template <GroupFamily F>
[[nodiscard]] CenteredOperations centered_operations(MatchContext<F> const &s,
                                                     Centering c) {
  CenteredOperations out;
  out.sym_translation.reserve(s.symmetry.size());
  out.sym_rotation.reserve(s.symmetry.size());
  out.db_translation.reserve(s.db_ops.size());
  for (auto const &op : s.symmetry) {
    out.sym_translation.push_back(transform_translation(c, op.translation));
    out.sym_rotation.push_back(transform_rotation(c, op.rotation));
  }
  for (auto const &op : s.db_ops) {
    out.db_translation.push_back(transform_translation(c, op.translation));
  }
  return out;
}

// Every operation must reproduce a distinct database operation once the
// origin shift is applied: same rotation, translations agreeing within
// symprec. Operations are matched in order, each to the first still-unmatched
// database operation carrying its rotation.
template <GroupFamily F>
[[nodiscard]] bool matches_database(MatchContext<F> const &s,
                                    CenteredOperations const &centered,
                                    Vector3d const &shift) {
  boost::container::static_vector<bool, 192> matched(s.db_ops.size(), false);
  for (auto const [i, op] : s.symmetry | std::views::enumerate) {
    auto const ui = static_cast<std::size_t>(i);
    Vector3d const lhs = centered.sym_translation[ui] + shift;
    Vector3d const shift_rot = centered.sym_rotation[ui] * shift;
    auto const [lo, hi] = s.db_by_rotation.equal_range(op.rotation);
    auto const it = std::ranges::find_if(lo, hi, [&](auto const &entry) {
      auto const idx = static_cast<std::size_t>(entry.second);
      return !matched[idx] &&
             coincident(lhs - centered.db_translation[idx], shift_rot,
                        s.primitive_lattice, s.symprec, s.periodicity());
    });
    if (it == hi) {
      return false;
    }
    matched[static_cast<std::size_t>(it->second)] = true;
  }
  return true;
}

// One (generators, VSpU) candidate: the origin shift if some choice of
// generator representatives makes the operations match the setting. The
// product varies the last generator fastest, so the first combination tried is
// the reference implementation's — the first operation carrying each rotation.
template <GroupFamily F>
[[nodiscard]] std::optional<Vector3d>
hall_symbol_shift(MatchContext<F> const &s, Centering c,
                  CenteredOperations const &centered,
                  data::GeneratorSet const &gens, data::VSpUSet const &vspu) {
  auto const choices = dw_choices(s, c, unpack_generators(gens));
  if (!choices) {
    return std::nullopt;
  }
  for (auto const &[dw0, dw1, dw2] : std::views::cartesian_product(
           (*choices)[0], (*choices)[1], (*choices)[2])) {
    Vector3d const shift = origin_shift(s, vspu, dw0, dw1, dw2);
    if (matches_database(s, centered, shift)) {
      return shift;
    }
  }
  return std::nullopt;
}

// The first candidate of a family that matches. The centering is fixed here,
// so this is where the shift-independent operation transforms are built.
template <GroupFamily F>
[[nodiscard]] std::optional<Vector3d>
first_shift(MatchContext<F> const &s, Centering c,
            std::span<data::GeneratorSet const> gens,
            std::span<data::VSpUSet const> vspu) {
  CenteredOperations const centered = centered_operations(s, c);
  for (auto const &[gen, vsp] : std::views::zip(gens, vspu)) {
    if (auto shift = hall_symbol_shift(s, c, centered, gen, vsp)) {
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

template <GroupFamily F>
[[nodiscard]] std::optional<Vector3d> dispatch(MatchContext<F> const &s,
                                               Centering c) {
  using namespace data;
  // Crystal system and the rhombohedral subsets are pure functions of the Hall
  // number, precomputed once in data::kHallClass. Layer groups (negative hall
  // numbers) carry no 3D hall-range bucket, so their crystal system is derived
  // from the point-group number; they are never rhombohedral and reuse the
  // same per-system generator families.
  HallClass const hc = [&] {
    if constexpr (F == GroupFamily::layer) {
      return HallClass{
          holohedry_from_pointgroup(spacegroup_type(s.hall).pointgroup_number),
          false, false};
    } else {
      return hall_class(s.hall.index());
    }
  }();

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

template <GroupFamily F>
std::optional<Vector3d>
SpacegroupMatcher<F>::match_hall(Matrix3d const &bravais_lattice,
                                 HallNumber hall, Centering centering,
                                 Operations const &symmetry, double symprec) {
  Operations const &db_ops = data::operations_from_database(hall);
  if (db_ops.size() != symmetry.size()) {
    return std::nullopt;
  }
  Matrix3d const primitive_lattice =
      bravais_lattice * centering_matrix_inv(centering);
  RotationMultimap<int> const symmetry_by_rotation =
      index_by_rotation(symmetry, &SymmetryOperation::rotation);
  MatchContext<F> const context{hall,     primitive_lattice,
                                symmetry, symmetry_by_rotation,
                                db_ops,   data::operations_by_rotation(hall),
                                symprec};
  // Bring the origin shift back to the bravais setting.
  return dispatch(context, centering).transform([&](Vector3d const &shift) {
    return Vector3d(centering_matrix_inv(centering) * shift);
  });
}

template std::optional<Vector3d>
SpacegroupMatcher<GroupFamily::space>::match_hall(Matrix3d const &, HallNumber,
                                                  Centering, Operations const &,
                                                  double);
template std::optional<Vector3d>
SpacegroupMatcher<GroupFamily::layer>::match_hall(Matrix3d const &, HallNumber,
                                                  Centering, Operations const &,
                                                  double);

} // namespace seitz::spacegroup
