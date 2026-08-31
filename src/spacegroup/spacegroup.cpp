#include <cppcrystal/spacegroup/spacegroup.hpp>

#include <cppcrystal/core/centering.hpp>
#include <cppcrystal/core/matrix_order.hpp>
#include <cppcrystal/core/point_group.hpp>
#include <cppcrystal/core/tolerance.hpp>
#include <cppcrystal/data/hall_classification.hpp>
#include <cppcrystal/math/fractional.hpp>
#include <cppcrystal/math/integer_matrix.hpp>
#include <cppcrystal/reduce/delaunay.hpp>
#include <cppcrystal/reduce/niggli.hpp>
#include <cppcrystal/spacegroup/hall_symbol.hpp>
#include <cppcrystal/symmetry/find_symmetry.hpp>
#include <cppcrystal/symmetry/pointgroup.hpp>

#include <algorithm>
#include <array>
#include <iterator>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

// Given a primitive cell's symmetry operations, find the Hall number,
// conventional/bravais lattice and origin shift: identify the point group, clean
// up the triclinic/monoclinic basis, determine the centering, build the
// conventional symmetry, then loop candidate Hall numbers permuting axis/setting
// choices and matching against the Hall database.
namespace cppcrystal::spacegroup {

using data::Centering;
using symmetry::Primitive;

namespace {

constexpr int kNumAttempt = 100;
constexpr double kReduceRate = 0.95;
constexpr double kIntPrec = 0.1;

// ---- matrix-of-doubles table builder ---------------------------------------

template <std::size_t N, std::size_t Rows = N>
[[nodiscard]] std::array<Matrix3d, N>
make_matrices(double const (&data)[Rows][3][3]) {
  std::array<Matrix3d, N> out;
  for (std::size_t i = 0; i < N; ++i) {
    out[i] << data[i][0][0], data[i][0][1], data[i][0][2], data[i][1][0],
        data[i][1][1], data[i][1][2], data[i][2][0], data[i][2][1],
        data[i][2][2];
  }
  return out;
}

// ---- change-of-basis tables ------------------------------------------------

// One axis/setting choice: the change-of-basis matrix, the centering that
// replaces C-centering under it, and the resulting unique axis. Stored as one
// row per choice so the three properties can never drift apart.
struct AxisChoice {
  Matrix3d basis;
  Centering centering;
  int unique_axis;
};

template <std::size_t N>
[[nodiscard]] std::array<AxisChoice, N>
make_axis_choices(double const (&bases)[N][3][3],
                  std::array<Centering, N> const &centerings,
                  std::array<int, N> const &unique_axes) {
  std::array<AxisChoice, N> out;
  auto const matrices = make_matrices<N>(bases);
  for (auto &&[choice, b, c, u] :
       std::views::zip(out, matrices, centerings, unique_axes)) {
    choice = {b, c, u};
  }
  return out;
}

// The two axes other than the unique one.
[[nodiscard]] auto non_unique_axes(int unique_axis) {
  return std::views::iota(0, 3) |
         std::views::filter([unique_axis](int j) { return j != unique_axis; });
}

// The 36 monoclinic change-of-basis choices (3D bulk; layer entries 36..47 are
// out of scope).
[[nodiscard]] std::array<AxisChoice, 36> const &monocli_axis_choices() {
  static double constexpr d[36][3][3] = {
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
      {{0, 0, 1}, {0, -1, 0}, {1, 0, 0}},
      {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
      {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
      {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
      {{-1, 0, 1}, {0, 1, 0}, {-1, 0, 0}},
      {{1, 0, -1}, {0, -1, 0}, {0, 0, -1}},
      {{0, 1, -1}, {1, 0, 0}, {0, 0, -1}},
      {{-1, -1, 0}, {0, 0, 1}, {-1, 0, 0}},
      {{1, -1, 0}, {0, 0, 1}, {0, -1, 0}},
      {{0, 1, 1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 0, -1}, {0, 1, 0}, {1, 0, -1}},
      {{-1, 0, 0}, {0, -1, 0}, {-1, 0, 1}},
      {{0, -1, 0}, {1, 0, 0}, {0, -1, 1}},
      {{0, 1, 0}, {0, 0, 1}, {1, 1, 0}},
      {{-1, 0, 0}, {0, 0, 1}, {-1, 1, 0}},
      {{0, 0, -1}, {1, 0, 0}, {0, -1, -1}},
      {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}},
      {{0, 0, -1}, {0, 1, 0}, {1, 0, 0}},
      {{0, 0, 1}, {-1, 0, 0}, {0, -1, 0}},
      {{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}},
      {{0, 1, 0}, {0, 0, -1}, {-1, 0, 0}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
      {{-1, 0, -1}, {0, -1, 0}, {-1, 0, 0}},
      {{1, 0, 1}, {0, 1, 0}, {0, 0, 1}},
      {{0, -1, -1}, {-1, 0, 0}, {0, 0, -1}},
      {{1, -1, 0}, {0, 0, -1}, {1, 0, 0}},
      {{-1, -1, 0}, {0, 0, -1}, {0, -1, 0}},
      {{0, -1, 1}, {-1, 0, 0}, {0, -1, 0}},
      {{0, 0, 1}, {0, -1, 0}, {1, 0, 1}},
      {{-1, 0, 0}, {0, 1, 0}, {-1, 0, -1}},
      {{0, 1, 0}, {-1, 0, 0}, {0, 1, 1}},
      {{0, 1, 0}, {0, 0, -1}, {-1, 1, 0}},
      {{1, 0, 0}, {0, 0, -1}, {1, 1, 0}},
      {{0, 0, -1}, {-1, 0, 0}, {0, 1, -1}},
  };
  using enum Centering;
  static constexpr std::array<Centering, 36> centerings = {
      c_face, a_face, b_face, b_face, a_face, c_face, a_face, c_face, c_face,
      a_face, b_face, b_face, body,   body,   body,   body,   body,   body,
      c_face, a_face, b_face, b_face, a_face, c_face, a_face, c_face, c_face,
      a_face, b_face, b_face, body,   body,   body,   body,   body,   body};
  static constexpr std::array<int, 36> unique_axes = {
      1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0,
      1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0, 1, 1, 0, 2, 2, 0};
  static auto const table = make_axis_choices<36>(d, centerings, unique_axes);
  return table;
}

// The 6 orthorhombic axis permutations {abc, cab, bca, ba-c, a-cb, -cba}.
[[nodiscard]] std::array<AxisChoice, 6> const &ortho_axis_choices() {
  static double constexpr d[6][3][3] = {
      {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},  {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
      {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},  {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},
      {{1, 0, 0}, {0, 0, 1}, {0, -1, 0}}, {{0, 0, 1}, {0, 1, 0}, {-1, 0, 0}}};
  using enum Centering;
  static constexpr std::array<Centering, 6> centerings = {
      c_face, b_face, a_face, c_face, b_face, a_face};
  static constexpr std::array<int, 6> unique_axes = {2, 1, 0, 2, 1, 0};
  static auto const table = make_axis_choices<6>(d, centerings, unique_axes);
  return table;
}

// Number of orthorhombic axis choices, indexed by (space-group number - 16),
// numbers 16..74.
[[nodiscard]] int num_axis_choices_ortho(int sg_number) {
  static constexpr std::array<int, 59> t = {
      6, 2, 2, 6, 2, 2, 6, 6, 6, 2, 1, 2, 1, 1, 1, 1, 2, 1, 2, 2,
      1, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 6, 6, 2, 2, 1, 1, 1, 1, 2,
      2, 1, 2, 2, 1, 3, 1, 1, 1, 2, 2, 2, 2, 6, 6, 6, 2, 6, 2};
  return t[static_cast<std::size_t>(sg_number - 16)];
}

// Number of orthorhombic layer-group axis choices, indexed by (layer-group
// number - 19); the 30 orthorhombic layer groups are numbers 19..48.
[[nodiscard]] int layer_num_axis_choices_ortho(int lg_number) {
  static constexpr std::array<int, 30> t = {
      2, 1,                         // 19-20
      2, 2, 2, 1, 2, 2, 1, 1, 1, 1, // 21-30
      1, 1, 1, 1, 1, 1, 2, 1, 2, 1, // 31-40
      1, 1, 1, 2, 1, 2, 2, 2};      // 41-48
  return t[static_cast<std::size_t>(lg_number - 19)];
}

// The number of orthorhombic axis choices for a setting, dispatching on the sign
// of the hall number (3D space group vs. layer group).
[[nodiscard]] int num_axis_choices_ortho_for(int hall_number, int number) {
  return hall_number < 0 ? layer_num_axis_choices_ortho(number)
                         : num_axis_choices_ortho(number);
}

[[nodiscard]] std::array<Matrix3d, 6> const &change_of_basis_rhombo() {
  static double const d[6][3][3] = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                    {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}},
                                    {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}},
                                    {{0, 0, -1}, {0, -1, 0}, {-1, 0, 0}},
                                    {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}},
                                    {{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}};
  static auto const table = make_matrices<6>(d);
  return table;
}

[[nodiscard]] std::array<Matrix3d, 6> const &change_of_basis_rhombo_hex() {
  static double const d[6][3][3] = {{{1, 0, 1}, {-1, 1, 1}, {0, -1, 1}},
                                    {{0, -1, 1}, {1, 0, 1}, {-1, 1, 1}},
                                    {{-1, 1, 1}, {0, -1, 1}, {1, 0, 1}},
                                    {{0, 1, -1}, {1, -1, -1}, {-1, 0, -1}},
                                    {{1, -1, -1}, {-1, 0, -1}, {0, 1, -1}},
                                    {{-1, 0, -1}, {0, 1, -1}, {1, -1, -1}}};
  static auto const table = make_matrices<6>(d);
  return table;
}

[[nodiscard]] std::array<Matrix3d, 8> const &change_of_basis_C4() {
  static double const d[8][3][3] = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                    {{0, -1, 0}, {1, 0, 0}, {0, 0, 1}},
                                    {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
                                    {{0, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
                                    {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},
                                    {{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}},
                                    {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}},
                                    {{1, 0, 0}, {0, -1, 0}, {0, 0, -1}}};
  static auto const table = make_matrices<8>(d);
  return table;
}

[[nodiscard]] std::array<Matrix3d, 12> const &change_of_basis_C6() {
  static double const d[12][3][3] = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                                     {{1, -1, 0}, {1, 0, 0}, {0, 0, 1}},
                                     {{0, -1, 0}, {1, -1, 0}, {0, 0, 1}},
                                     {{-1, 0, 0}, {0, -1, 0}, {0, 0, 1}},
                                     {{-1, 1, 0}, {-1, 0, 0}, {0, 0, 1}},
                                     {{0, 1, 0}, {-1, 1, 0}, {0, 0, 1}},
                                     {{0, 1, 0}, {1, 0, 0}, {0, 0, -1}},
                                     {{-1, 1, 0}, {0, 1, 0}, {0, 0, -1}},
                                     {{-1, 0, 0}, {-1, 1, 0}, {0, 0, -1}},
                                     {{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}},
                                     {{1, -1, 0}, {0, -1, 0}, {0, 0, -1}},
                                     {{1, 0, 0}, {1, -1, 0}, {0, 0, -1}}};
  static auto const table = make_matrices<12>(d);
  return table;
}

// ---- correction / fixed matrices -------------------------------------------

[[nodiscard]] Matrix3d const &monocli_i2c() {
  static Matrix3d const m =
      (Matrix3d() << 1, 0, -1, 0, 1, 0, 1, 0, 0).finished();
  return m;
}
[[nodiscard]] Matrix3d const &monocli_a2c() {
  static Matrix3d const m =
      (Matrix3d() << 0, 0, 1, 0, -1, 0, 1, 0, 0).finished();
  return m;
}
[[nodiscard]] Matrix3d const &a2c() {
  static Matrix3d const m =
      (Matrix3d() << 0, 0, 1, 1, 0, 0, 0, 1, 0).finished();
  return m;
}
[[nodiscard]] Matrix3d const &b2c() {
  static Matrix3d const m =
      (Matrix3d() << 0, 1, 0, 0, 0, 1, 1, 0, 0).finished();
  return m;
}
[[nodiscard]] Matrix3d const &rhombo_obverse() {
  static Matrix3d const m = (Matrix3d() << 2. / 3, -1. / 3, -1. / 3, 1. / 3,
                             1. / 3, -2. / 3, 1. / 3, 1. / 3, 1. / 3)
                                .finished();
  return m;
}
[[nodiscard]] Matrix3d const &rhomb_reverse() {
  static Matrix3d const m = (Matrix3d() << 1. / 3, -2. / 3, 1. / 3, 2. / 3,
                             -1. / 3, -1. / 3, 1. / 3, 1. / 3, 1. / 3)
                                .finished();
  return m;
}

// spacegroup_to_hall_number[230]: first Hall number of each space group.
[[nodiscard]] std::span<int const> spacegroup_to_hall_number() {
  static constexpr std::array<int, 230> t = {
      1,   2,   3,   6,   9,   18,  21,  30,  39,  57,  60,  63,  72,  81,  90,
      108, 109, 112, 115, 116, 119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
      155, 161, 164, 170, 173, 176, 182, 185, 191, 197, 203, 209, 212, 215, 218,
      221, 227, 228, 230, 233, 239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
      290, 292, 298, 304, 310, 313, 316, 322, 334, 335, 337, 338, 341, 343, 349,
      350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
      368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381, 382,
      383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
      398, 399, 400, 401, 402, 404, 406, 407, 408, 410, 412, 413, 414, 416, 418,
      419, 420, 422, 424, 425, 426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
      440, 441, 442, 443, 444, 446, 447, 448, 449, 450, 452, 454, 455, 456, 457,
      458, 460, 462, 463, 464, 465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
      475, 476, 477, 478, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489,
      490, 491, 492, 493, 494, 495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
      507, 508, 509, 510, 511, 512, 513, 514, 515, 516, 517, 518, 520, 521, 523,
      524, 525, 527, 529, 530};
  return t;
}

// layer_group_to_hall_number[116]: the (negative) layer Hall settings searched
// for a layer cell, one representative per setting.
[[nodiscard]] std::span<int const> layer_group_to_hall_number() {
  static constexpr std::array<int, 80> t = {
      -1,   -2,   -3,   -4,   -5,   -8,   -9,   -12,  -14,  -16,  -18,  -20,
      -22,  -24,  -26,  -28,  -30,  -32,  -34,  -35,  -37,  -38,  -39,  -40,
      -42,  -43,  -44,  -46,  -48,  -50,  -52,  -54,  -56,  -58,  -60,  -62,
      -64,  -65,  -67,  -68,  -70,  -72,  -74,  -76,  -77,  -79,  -80,  -81,
      -82,  -83,  -84,  -85,  -87,  -88,  -89,  -90,  -91,  -92,  -93,  -94,
      -95,  -96,  -98,  -99,  -101, -102, -103, -104, -105, -106, -107, -108,
      -109, -110, -111, -112, -113, -114, -115, -116};
  return t;
}

// ---- small matrix predicates -----------------------------------------------

// |a - b| <= symprec elementwise.
[[nodiscard]] bool matrices_close(Matrix3d const &a, Matrix3d const &b,
                                  double symprec) {
  return (a - b).cwiseAbs().maxCoeff() <= symprec;
}

// ---- is_equivalent_lattice -------------------------------------------------

// Isometric transform tmat such that orig = lattice * tmat, if the two lattices
// are equivalent under `mode`; nullopt otherwise.
//   mode 0: tmat == identity
//   mode 1: |tmat| == identity (axis flips allowed)
//   mode 2: tmat integer & unimodular & metric tensors agree
[[nodiscard]] std::optional<Matrix3d>
is_equivalent_lattice(int mode, Matrix3d const &lattice, Matrix3d const &orig,
                      double symprec) {
  if (std::abs(lattice.determinant() - orig.determinant()) > symprec) {
    return std::nullopt;
  }
  auto const inv = math::inverse(lattice, symprec);
  if (!inv) {
    return std::nullopt;
  }
  Matrix3d const tmat = *inv * orig;

  switch (mode) {
  case 0:
    if (matrices_close(Matrix3d::Identity(), tmat, symprec)) {
      return tmat;
    }
    break;
  case 1:
    if (matrices_close(Matrix3d::Identity(), tmat.cwiseAbs(), symprec)) {
      return tmat;
    }
    break;
  case 2: {
    Matrix3i const tmat_int = math::round_to_int(tmat);
    if (!matrices_close(tmat_int.cast<double>(), tmat, symprec)) {
      break;
    }
    if (tmat_int.determinant() != 1) {
      break;
    }
    if (matrices_close(math::metric_tensor(orig), math::metric_tensor(lattice),
                       symprec)) {
      return tmat;
    }
    break;
  }
  default:
    break;
  }
  return std::nullopt;
}

// ---- centering -------------------------------------------------------------

// Detect the base-centering type encoded in an integer transformation matrix.
[[nodiscard]] Centering get_base_center(Matrix3i const &tmat) {
  auto const axis = std::views::iota(0, 3);
  if (std::ranges::any_of(axis, [&](int i) { // C center
        return tmat(i, 0) == 0 && tmat(i, 1) == 0 && std::abs(tmat(i, 2)) == 1;
      })) {
    return Centering::c_face;
  }
  if (std::ranges::any_of(axis, [&](int i) { // A center
        return std::abs(tmat(i, 0)) == 1 && tmat(i, 1) == 0 && tmat(i, 2) == 0;
      })) {
    return Centering::a_face;
  }
  if (std::ranges::any_of(axis, [&](int i) { // B center
        return tmat(i, 0) == 0 && std::abs(tmat(i, 1)) == 1 && tmat(i, 2) == 0;
      })) {
    return Centering::b_face;
  }
  // body center: every row's absolute coordinate sum is 2
  bool const body = std::ranges::all_of(axis, [&](int i) {
    return std::abs(tmat(i, 0)) + std::abs(tmat(i, 1)) + std::abs(tmat(i, 2)) ==
           2;
  });
  return body ? Centering::body : Centering::primitive;
}

struct CenteringResult {
  Centering centering;
  Matrix3d correction_mat;
};

// Centering type + correction matrix from the determinant (multiplicity) of the
// integer transformation matrix and the Laue class.
[[nodiscard]] std::optional<CenteringResult> get_centering(Matrix3i const &tmat,
                                                           Laue laue) {
  Matrix3d correction = Matrix3d::Identity();
  int const det = std::abs(tmat.determinant());

  switch (det) {
  case 1:
    return CenteringResult{Centering::primitive, correction};
  case 2: {
    Centering centering = get_base_center(tmat);
    if (centering == Centering::a_face) {
      correction = (laue == Laue::laue_2m) ? monocli_a2c() : a2c();
      centering = Centering::c_face;
    } else if (centering == Centering::b_face) {
      correction = b2c();
      centering = Centering::c_face;
    } else if (laue == Laue::laue_2m && centering == Centering::body) {
      correction = monocli_i2c();
      centering = Centering::c_face;
    }
    return CenteringResult{centering, correction};
  }
  case 3: { // hP (a=b) but not hR (a=b=c)
    if (math::is_int_matrix(Matrix3d(tmat.cast<double>() * rhombo_obverse()),
                            kIntPrec)) {
      correction = rhombo_obverse();
    }
    if (math::is_int_matrix(Matrix3d(tmat.cast<double>() * rhomb_reverse()),
                            kIntPrec)) {
      correction = rhomb_reverse();
    }
    return CenteringResult{Centering::r_center, correction};
  }
  case 4:
    return CenteringResult{Centering::face, correction};
  default:
    return std::nullopt;
  }
}

// centering_shifts (the non-trivial centering translations) lives in
// core/centering.hpp, shared with cell standardization.

// ---- conventional symmetry -------------------------------------------------

// Recover the conventional-cell operations from the primitive ones via the
// similarity transform C S C^-1 (C = tmat), then replicate them across the
// centering translations.
[[nodiscard]] SymmetryOperations
get_conventional_symmetry(Matrix3d const &tmat, Centering centering,
                          SymmetryOperations const &primitive_sym) {
  Matrix3d const inv_tmat = tmat.inverse();

  SymmetryOperations base;
  base.reserve(primitive_sym.size());
  std::ranges::transform(primitive_sym, std::back_inserter(base),
                         [&](SymmetryOperation const &op) {
                           return conjugated_by(op, inv_tmat, tmat);
                         });

  SymmetryOperations out = base;
  for (Vector3d const &shift : centering_shifts(centering)) {
    for (auto const &op : base) {
      out.push_back({op.rotation, Vector3d(op.translation + shift)});
    }
  }
  return out;
}

// Rhombohedral (R_CENTER) keeps the primitive a=b=c basis, so its conventional
// symmetry is built as PRIMITIVE.
[[nodiscard]] SymmetryOperations
get_initial_conventional_symmetry(Centering centering, Matrix3d const &tmat,
                                  SymmetryOperations const &symmetry) {
  Centering const c =
      (centering == Centering::r_center) ? Centering::primitive : centering;
  return get_conventional_symmetry(tmat, c, symmetry);
}

// ---- triclinic / monoclinic basis cleanup ----------------------------------

// Niggli-reduce the conventional lattice and update the integer transformation
// matrix to the smallest right-handed basis.
[[nodiscard]] std::optional<Matrix3i>
change_basis_tricli(Matrix3d const &conv_lattice,
                    Matrix3d const &primitive_lattice, double symprec) {
  auto const reduced = reduce::niggli_reduce(conv_lattice, symprec * symprec);
  if (!reduced) {
    return std::nullopt;
  }
  Matrix3d const smallest =
      reduced->determinant() < 0 ? Matrix3d(-*reduced) : *reduced;
  Matrix3d const tmat = primitive_lattice.inverse() * smallest;
  return math::round_to_int(tmat);
}

// 2D-Delaunay-reduce the plane orthogonal to the kept axis and update the
// integer transformation matrix. For a 3D monoclinic cell the kept axis is the
// unique axis b; for a layer the kept axis is the aperiodic axis (in its
// conventional position), so the reduction acts on the periodic plane.
[[nodiscard]] std::optional<Matrix3i>
change_basis_monocli(Matrix3d const &conv_lattice,
                     Matrix3d const &primitive_lattice, double symprec,
                     std::optional<int> aperiodic_axis_conv) {
  int const keep_axis = aperiodic_axis_conv.value_or(1);
  auto const smallest =
      reduce::delaunay_reduce(conv_lattice, keep_axis, symprec);
  if (!smallest) {
    return std::nullopt;
  }
  Matrix3d const tmat = primitive_lattice.inverse() * *smallest;
  return math::round_to_int(tmat);
}

// ---- match_hall_symbol_db family -------------------------------------------

// A successful per-system match: the (updated) conventional lattice and the
// origin shift aligning the operations with the Hall database.
struct MatchResult {
  Vector3d origin_shift;
  Matrix3d conv_lattice;
};

// Thin wrapper over match_hall_symbol returning a MatchResult.
[[nodiscard]] std::optional<MatchResult>
try_hall(Matrix3d const &conv_lattice, int hall_number, Centering centering,
         SymmetryOperations const &symmetry, double symprec) {
  auto const shift = match_hall_symbol(conv_lattice, hall_number, centering,
                                       symmetry, symprec);
  if (!shift) {
    return std::nullopt;
  }
  return MatchResult{*shift, conv_lattice};
}

// The two-pass preference loop shared by the per-system matchers: when the
// input lattice is usable, first try every choice constrained to it (so a
// basis matching the input wins), then retry unconstrained.
// `attempt(choice, orig_lattice_or_null)` produces the candidate result.
template <std::ranges::input_range R, class Try>
[[nodiscard]] std::optional<MatchResult>
first_match(R &&choices, Matrix3d const *orig_lattice, double symprec,
            Try &&attempt) {
  if (orig_lattice && orig_lattice->determinant() > symprec) {
    for (auto const &choice : choices) {
      if (auto r = attempt(choice, orig_lattice)) {
        return r;
      }
    }
  }
  for (auto const &choice : choices) {
    if (auto r = attempt(choice, nullptr)) {
      return r;
    }
  }
  return std::nullopt;
}

// Shared loop for TETRA / HEXA / TRIGO / rhombohedral: try each rotation choice
// (with the input-similarity preference first), defer to the Hall matcher.
[[nodiscard]] std::optional<MatchResult> match_db_change_of_basis_loop(
    Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
    std::span<Matrix3d const> change_of_basis, int hall_number,
    Centering centering, SymmetryOperations const &conv_symmetry,
    double symprec) {
  Centering const centering_for_symmetry = (centering == Centering::r_center)
                                               ? Centering::r_center
                                               : Centering::primitive;

  return first_match(
      change_of_basis, orig_lattice, symprec,
      [&](Matrix3d const &cob,
          Matrix3d const *orig) -> std::optional<MatchResult> {
        Matrix3d const changed_lattice = conv_lattice * cob;
        if (orig &&
            !is_equivalent_lattice(0, changed_lattice, *orig, symprec)) {
          return std::nullopt;
        }
        SymmetryOperations const changed_symmetry = get_conventional_symmetry(
            cob, centering_for_symmetry, conv_symmetry);
        return try_hall(changed_lattice, hall_number, centering,
                        changed_symmetry, symprec);
      });
}

[[nodiscard]] std::optional<MatchResult>
match_db_rhombo(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
                int hall_number, SymmetryOperations const &conv_symmetry,
                double symprec) {
  if (data::is_rhombo_hex_setting(hall_number)) { // hexagonal (hP) setting
    return match_db_change_of_basis_loop(
        conv_lattice, orig_lattice, change_of_basis_rhombo_hex(), hall_number,
        Centering::r_center, conv_symmetry, symprec);
  }
  // rhombohedral (a=b=c) setting
  return match_db_change_of_basis_loop(
      conv_lattice, orig_lattice, change_of_basis_rhombo(), hall_number,
      Centering::primitive, conv_symmetry, symprec);
}

[[nodiscard]] std::optional<MatchResult>
match_db_others(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
                int hall_number, Centering centering, Holohedry holohedry,
                SymmetryOperations const &conv_symmetry, double symprec) {
  if (holohedry == Holohedry::triclinic) {
    return try_hall(conv_lattice, hall_number, centering, conv_symmetry,
                    symprec);
  }
  if (holohedry == Holohedry::tetragonal) {
    return match_db_change_of_basis_loop(conv_lattice, orig_lattice,
                                         change_of_basis_C4(), hall_number,
                                         centering, conv_symmetry, symprec);
  }
  // hexagonal / trigonal (non-rhombohedral)
  return match_db_change_of_basis_loop(conv_lattice, orig_lattice,
                                       change_of_basis_C6(), hall_number,
                                       centering, conv_symmetry, symprec);
}

[[nodiscard]] std::optional<MatchResult> match_db_cubic_in_loop(
    Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
    Matrix3d const &change, int hall_number, Centering centering,
    SymmetryOperations const &conv_symmetry, double symprec) {
  Matrix3d change_of_basis = change;
  Matrix3d changed_lattice = conv_lattice * change_of_basis;

  if (orig_lattice) {
    auto const tmat =
        is_equivalent_lattice(1, changed_lattice, *orig_lattice, symprec);
    if (!tmat) {
      return std::nullopt;
    }
    changed_lattice = changed_lattice * *tmat;
    change_of_basis = change_of_basis * *tmat;
  }

  SymmetryOperations const changed_symmetry = get_conventional_symmetry(
      change_of_basis, Centering::primitive, conv_symmetry);
  return try_hall(changed_lattice, hall_number, centering, changed_symmetry,
                  symprec);
}

[[nodiscard]] std::optional<MatchResult>
match_db_cubic(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
               int hall_number, Centering centering,
               SymmetryOperations const &conv_symmetry, double symprec) {
  return first_match(ortho_axis_choices(), orig_lattice, symprec,
                     [&](AxisChoice const &choice, Matrix3d const *orig) {
                       return match_db_cubic_in_loop(
                           conv_lattice, orig, choice.basis, hall_number,
                           centering, conv_symmetry, symprec);
                     });
}

// Norm-ordering preference among the free axes. Returns false if the choice
// violates the required |a|<=|b|<=|c| order.
[[nodiscard]] bool ortho_axis_norms_ok(Matrix3d const &lattice, int unique_axis,
                                       int num_free_axes) {
  auto sqnorm = [&](int j) { return lattice.col(j).squaredNorm(); };
  if (num_free_axes == 2) {
    auto norms = non_unique_axes(unique_axis) | std::views::transform(sqnorm);
    auto it = norms.begin();
    double const first = *it;
    return !sqnorm_longer(first, *++it);
  }
  if (num_free_axes == 3) {
    return !(sqnorm_longer(sqnorm(0), sqnorm(1)) ||
             sqnorm_longer(sqnorm(0), sqnorm(2)));
  }
  if (num_free_axes == 6) {
    return !(sqnorm_longer(sqnorm(0), sqnorm(1)) ||
             sqnorm_longer(sqnorm(1), sqnorm(2)));
  }
  return true; // num_free_axes == 0 (principal-axis search) or 1
}

[[nodiscard]] std::optional<MatchResult> match_db_ortho_in_loop(
    Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
    AxisChoice const &choice, int hall_number, Centering centering,
    SymmetryOperations const &symmetry, int num_free_axes, double symprec) {
  Centering const changed_centering =
      (centering == Centering::c_face) ? choice.centering : centering;

  Matrix3d change_of_basis = choice.basis;
  Matrix3d changed_lattice = conv_lattice * change_of_basis;

  if (orig_lattice) {
    auto const tmat =
        is_equivalent_lattice(1, changed_lattice, *orig_lattice, symprec);
    if (!tmat) {
      return std::nullopt;
    }
    changed_lattice = changed_lattice * *tmat;
    change_of_basis = change_of_basis * *tmat;
  }

  if (!ortho_axis_norms_ok(changed_lattice, choice.unique_axis,
                           num_free_axes)) {
    return std::nullopt;
  }

  SymmetryOperations const changed_symmetry = get_conventional_symmetry(
      change_of_basis, Centering::primitive, symmetry);
  return try_hall(changed_lattice, hall_number, changed_centering,
                  changed_symmetry, symprec);
}

[[nodiscard]] std::optional<MatchResult>
match_db_ortho(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
               int hall_number, Centering centering,
               SymmetryOperations const &symmetry, int num_free_axes,
               double symprec) {
  // 3D tries all six axis permutations {abc, bca, cab, ba-c, a-cb, -cba}; a
  // layer group (hall < 0) tries only abc and ba-c (stride 3), the two that
  // keep the aperiodic axis at c.
  int const step = hall_number < 0 ? 3 : 1;
  auto const choices = std::views::iota(0, 6) | std::views::stride(step) |
                       std::views::transform([](int i) -> AxisChoice const & {
                         return ortho_axis_choices()[
                             static_cast<std::size_t>(i)];
                       });
  return first_match(choices, orig_lattice, symprec,
                     [&](AxisChoice const &choice, Matrix3d const *orig) {
                       return match_db_ortho_in_loop(
                           conv_lattice, orig, choice, hall_number, centering,
                           symmetry, num_free_axes, symprec);
                     });
}

// One change-of-basis attempt for monoclinic. Status: 0 not found, 1 found, 2
// found and the basis matches the input lattice's (a,b,c) choice.
struct MonocliCandidate {
  int status = 0;
  Vector3d origin_shift{Vector3d::Zero()};
  Matrix3d conv_lattice{Matrix3d::Identity()};
  std::array<double, 2> norms_squared{0.0, 0.0};
};

[[nodiscard]] MonocliCandidate
match_db_monocli_in_loop(Matrix3d conv_lattice, AxisChoice const &choice,
                         Matrix3d const *orig_lattice, bool check_norms,
                         int hall_number, Centering centering,
                         SymmetryOperations const &conv_symmetry,
                         double symprec) {
  Centering const changed_centering =
      (centering == Centering::c_face) ? choice.centering : centering;

  Matrix3d change_of_basis = choice.basis;
  conv_lattice = conv_lattice * change_of_basis;
  int const unique_axis = choice.unique_axis;

  // The two non-unique axes and their squared norms.
  std::array<Vector3d, 2> vec;
  std::array<double, 2> norms_squared;
  for (int l = 0; int const j : non_unique_axes(unique_axis)) {
    vec[static_cast<std::size_t>(l)] = conv_lattice.col(j);
    norms_squared[static_cast<std::size_t>(l)] = vec[static_cast<std::size_t>(l)].squaredNorm();
    ++l;
  }

  // Discard if the principal angle is acute.
  if (vec[0].dot(vec[1]) > kZeroPrec) {
    return {};
  }
  // Prefer |a| <= |b| <= |c| among the non-principal axes when free to choose.
  if (check_norms && sqnorm_longer(norms_squared[0], norms_squared[1])) {
    return {};
  }

  int status = 1;
  if (orig_lattice && orig_lattice->determinant() > symprec) {
    // mode-1 equivalence effectively checks the C2 rotation about the unique
    // axis; only accept flips that keep the unique axis fixed.
    if (auto const tmat =
            is_equivalent_lattice(1, conv_lattice, *orig_lattice, symprec)) {
      int const u1 = (unique_axis + 1) % 3;
      int const u2 = (unique_axis + 2) % 3;
      if ((*tmat)(u1, u1) * (*tmat)(u2, u2) > kZeroPrec) {
        conv_lattice = conv_lattice * *tmat;
        change_of_basis = change_of_basis * *tmat;
        status = 2;
      }
    }
  }

  SymmetryOperations const changed_symmetry = get_conventional_symmetry(
      change_of_basis, Centering::primitive, conv_symmetry);
  auto const shift = match_hall_symbol(
      conv_lattice, hall_number, changed_centering, changed_symmetry, symprec);
  if (!shift) {
    return {};
  }
  return {status, *shift, conv_lattice, norms_squared};
}

[[nodiscard]] std::optional<MatchResult>
match_db_monocli(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
                 int hall_number, int group_number, Centering centering,
                 SymmetryOperations const &conv_symmetry, double symprec) {
  // E. Parthe & L. M. Gelato (1983): the cell preference is enforced (a/b/c
  // length ordering) only for the 5 monoclinic space-group types below.
  bool const check_norms = group_number == 3 || group_number == 4 ||
                           group_number == 6 || group_number == 10 ||
                           group_number == 11;

  std::array<MonocliCandidate, 36> found;
  std::ranges::transform(monocli_axis_choices(), found.begin(),
                         [&](AxisChoice const &choice) {
                           return match_db_monocli_in_loop(
                               conv_lattice, choice, orig_lattice, check_norms,
                               hall_number, centering, conv_symmetry, symprec);
                         });

  auto is_found = [](MonocliCandidate const &c) { return c.status != 0; };
  if (std::ranges::none_of(found, is_found)) {
    return std::nullopt;
  }

  auto norm_sum = [](MonocliCandidate const &c) {
    return std::sqrt(c.norms_squared[0]) + std::sqrt(c.norms_squared[1]);
  };

  // Shortest pair of non-unique axes among the matches (at least one exists,
  // per the none_of guard above).
  double const shortest = std::ranges::min(
      found | std::views::filter(is_found) | std::views::transform(norm_sum));

  // Among the near-shortest matches, prefer a basis that matches the input
  // (status == 2); max_element keeps the first of equal status, matching the
  // historical first-candidate preference.
  auto near_shortest =
      found | std::views::filter(is_found) |
      std::views::filter([&](MonocliCandidate const &c) {
        return std::abs(norm_sum(c) - shortest) < symprec;
      });
  auto const chosen =
      std::ranges::max_element(near_shortest, {}, &MonocliCandidate::status);
  return MatchResult{chosen->origin_shift, chosen->conv_lattice};
}

// Dispatch by holohedry to the per-system matcher.
[[nodiscard]] std::optional<MatchResult>
match_hall_symbol_db(Matrix3d const &conv_lattice, Matrix3d const *orig_lattice,
                     int hall_number, int pointgroup_number,
                     Holohedry holohedry, Centering centering,
                     SymmetryOperations const &symmetry, double symprec) {
  data::SpacegroupType const sg = data::spacegroup_type(hall_number);
  if (pointgroup_number != sg.pointgroup_number) {
    return std::nullopt;
  }

  switch (holohedry) {
  case Holohedry::monoclinic:
    return match_db_monocli(conv_lattice, orig_lattice, hall_number, sg.number,
                            centering, symmetry, symprec);

  case Holohedry::orthorhombic: {
    int const num_free_axes = num_axis_choices_ortho_for(hall_number, sg.number);
    // Two-axis case: first fix the principal axis for the representative Hall
    // symbol, then match the requested one in the resulting basis. This
    // representative-Hall two-step is a 3D-only refinement; layer groups match
    // directly with the c-preserving axis choices.
    if (hall_number > 0 && num_free_axes == 2) {
      auto const rep = match_db_ortho(
          conv_lattice, orig_lattice,
          spacegroup_to_hall_number()[static_cast<std::size_t>(sg.number - 1)],
          centering, symmetry, 0, symprec);
      if (!rep) {
        return std::nullopt;
      }
      Matrix3d const changed_lattice = rep->conv_lattice;
      Matrix3d const tmat = conv_lattice.inverse() * changed_lattice;
      SymmetryOperations const changed_symmetry =
          get_conventional_symmetry(tmat, Centering::primitive, symmetry);
      return match_db_ortho(changed_lattice, orig_lattice, hall_number,
                            centering, changed_symmetry, 2, symprec);
    }
    return match_db_ortho(conv_lattice, orig_lattice, hall_number, centering,
                          symmetry, num_free_axes, symprec);
  }

  case Holohedry::cubic:
    return match_db_cubic(conv_lattice, orig_lattice, hall_number, centering,
                          symmetry, symprec);

  case Holohedry::trigonal:
    // Rhombohedral subset (R-centered) uses the a=b=c / hP machinery.
    if (centering == Centering::r_center &&
        data::is_rhombohedral_hall(hall_number)) {
      return match_db_rhombo(conv_lattice, orig_lattice, hall_number, symmetry,
                             symprec);
    }
    [[fallthrough]]; // other trigonal cases: shared change-of-basis loop
  default:           // hexagonal, tetragonal, triclinic, rest of trigonal
    return match_db_others(conv_lattice, orig_lattice, hall_number, centering,
                           holohedry, symmetry, symprec);
  }
}

// ---- search ----------------------------------------------------------------

struct SearchResult {
  int hall_number;
  Matrix3d conv_lattice;
  Vector3d origin_shift;
};

// For the given operations, find the conventional setting and the first
// matching Hall number among the candidates.
[[nodiscard]] std::optional<SearchResult>
search_hall_number(std::span<int const> candidates, Primitive const &primitive,
                   SymmetryOperations const &symmetry, double symprec) {
  std::vector<Matrix3i> const rotations = symmetry::rotations_of(symmetry);
  auto const ptg =
      symmetry::get_pointgroup(rotations, primitive.cell.aperiodic_axis());
  if (!ptg || ptg->pointgroup.number == 0) {
    return std::nullopt;
  }

  Matrix3d const &prim_lat = primitive.cell.lattice();
  Matrix3i tmat_int = ptg->transformation;
  Laue const laue = ptg->pointgroup.laue;

  // For LAUE1 / LAUE2M, make the smallest lattice (tricli: Niggli; monocli: 2D
  // Delaunay) and update the integer transformation.
  if (laue == Laue::laue_1 || laue == Laue::laue_2m) {
    Matrix3d const conv_tmp = prim_lat * tmat_int.cast<double>();
    // For a layer, locate the aperiodic axis in the conventional setting (the
    // column the primitive aperiodic axis maps to under tmat_int) so the
    // monoclinic reduction keeps it and reduces the periodic plane.
    std::optional<int> aperiodic_conv;
    if (auto const ap = primitive.cell.aperiodic_axis()) {
      // The last nonzero column wins, matching the historical scan direction.
      auto const cols = std::views::iota(0, 3) | std::views::reverse;
      auto const it = std::ranges::find_if(
          cols, [&](int i) { return tmat_int(*ap, i) != 0; });
      if (it != cols.end()) {
        aperiodic_conv = *it;
      }
    }
    auto const updated =
        (laue == Laue::laue_1)
            ? change_basis_tricli(conv_tmp, prim_lat, symprec)
            : change_basis_monocli(conv_tmp, prim_lat, symprec, aperiodic_conv);
    if (!updated) {
      return std::nullopt;
    }
    tmat_int = *updated;
  }

  auto const centering = get_centering(tmat_int, laue);
  if (!centering) {
    return std::nullopt;
  }

  Matrix3d const tmat = tmat_int.cast<double>() * centering->correction_mat;
  Matrix3d const conv_lattice = prim_lat * tmat;

  SymmetryOperations const conv_symmetry =
      get_initial_conventional_symmetry(centering->centering, tmat, symmetry);

  Matrix3d const &orig_lattice = primitive.orig_lattice;
  for (int const hall_number : candidates) {
    auto const match =
        match_hall_symbol_db(conv_lattice, &orig_lattice, hall_number,
                             ptg->pointgroup.number, ptg->pointgroup.holohedry,
                             centering->centering, conv_symmetry, symprec);
    if (match) {
      return SearchResult{hall_number, match->conv_lattice,
                          match->origin_shift};
    }
  }
  return std::nullopt;
}

// Try once, then progressively tighten the operation set (reduce_symmetry)
// until a Hall number is found.
[[nodiscard]] std::optional<SearchResult>
iterative_search_hall_number(std::span<int const> candidates,
                             Primitive const &primitive,
                             SymmetryOperations const &symmetry, double symprec,
                             AngleTolerance angle_tolerance) {
  if (auto r = search_hall_number(candidates, primitive, symmetry, symprec)) {
    return r;
  }

  double tolerance = symprec;
  for (int attempt = 0; attempt < kNumAttempt; ++attempt) {
    tolerance *= kReduceRate;
    SymmetryOperations const reduced = symmetry::reduce_symmetry(
        primitive.cell, symmetry, tolerance, angle_tolerance);
    if (reduced.empty()) {
      continue;
    }
    if (auto r = search_hall_number(candidates, primitive, reduced, symprec)) {
      return r;
    }
  }
  return std::nullopt;
}

// Sanity check: in a primitive cell every operation must carry a distinct
// rotation.
[[nodiscard]] bool point_symmetry_intact(SymmetryOperations const &symmetry) {
  return unique_by_rotation(symmetry, &SymmetryOperation::rotation).size() ==
         symmetry.size();
}

[[nodiscard]] Result<Spacegroup> search_spacegroup_with_symmetry(
    std::span<int const> candidates, Primitive const &primitive,
    SymmetryOperations const &symmetry, double symprec,
    AngleTolerance angle_tolerance) {
  if (!point_symmetry_intact(symmetry)) {
    return leaf::new_error(e_spacegroup_search_failed{});
  }

  auto const found = iterative_search_hall_number(
      candidates, primitive, symmetry, symprec, angle_tolerance);
  if (!found) {
    return leaf::new_error(e_spacegroup_search_failed{});
  }

  return Spacegroup{data::spacegroup_type(found->hall_number),
                    found->conv_lattice, found->origin_shift};
}

} // namespace

Result<Spacegroup> search_spacegroup(Primitive const &primitive,
                                     int hall_number, double symprec,
                                     AngleTolerance angle_tolerance) {
  BOOST_LEAF_AUTO(symmetry, symmetry::find_symmetry(primitive.cell, symprec,
                                                    angle_tolerance));

  if (hall_number != 0) {
    std::array<int, 1> const candidate = {hall_number};
    return search_spacegroup_with_symmetry(candidate, primitive, symmetry,
                                           symprec, angle_tolerance);
  }
  // A layer cell (aperiodic axis set) searches the 80 layer-group settings;
  // a 3D cell searches the 230 space groups.
  std::span<int const> const candidates =
      primitive.cell.aperiodic_axis() ? layer_group_to_hall_number()
                                      : spacegroup_to_hall_number();
  return search_spacegroup_with_symmetry(candidates, primitive, symmetry,
                                         symprec, angle_tolerance);
}

Result<Spacegroup>
search_spacegroup_with_symmetry(SymmetryOperations const &operations,
                                Matrix3d const &prim_lattice, double symprec) {
  // A single notional atom at the origin; only the lattice matters here.
  Positions pos(1, 3);
  pos.setZero();
  symmetry::Primitive const primitive{Cell(prim_lattice, pos, Types{1}),
                                      std::vector<int>{0}, prim_lattice,
                                      symprec, std::nullopt};
  return search_spacegroup_with_symmetry(spacegroup_to_hall_number(), primitive,
                                         operations, symprec, std::nullopt);
}

template <LatticeSetting Setting>
Result<Spacegroup>
spacegroup_type_from_symmetry(SymmetryOperations const &operations,
                              Matrix3d const &lattice, double symprec) {
  auto const prim = symmetry::primitive_symmetry(operations, symprec);
  if (!prim) {
    return leaf::new_error(e_spacegroup_search_failed{});
  }
  SymmetryOperations const &prim_sym = prim->first;
  Matrix3d const &t_mat = prim->second;

  Matrix3d prim_lat;
  if constexpr (Setting == LatticeSetting::conventional) {
    prim_lat = lattice * t_mat.inverse();
  } else {
    prim_lat = lattice;
  }

  // Niggli-reduce the primitive cell (required by the Hall-symbol matcher),
  // then bring the operations into that reduced basis (the rotations are already
  // distinct, so it is just the change of basis).
  BOOST_LEAF_AUTO(red_lat, reduce::niggli_reduce(prim_lat, symprec));
  Matrix3d const t_mat2 = red_lat.inverse() * prim_lat;
  Matrix3d const inv2 = t_mat2.inverse();
  SymmetryOperations red_sym;
  red_sym.reserve(prim_sym.size());
  std::ranges::transform(prim_sym, std::back_inserter(red_sym),
                         [&](SymmetryOperation const &op) {
                           return conjugated_by(op, t_mat2, inv2);
                         });

  return search_spacegroup_with_symmetry(red_sym, red_lat, symprec);
}

// Both lattice-setting specializations are used.
template Result<Spacegroup>
spacegroup_type_from_symmetry<LatticeSetting::conventional>(
    SymmetryOperations const &, Matrix3d const &, double);
template Result<Spacegroup>
spacegroup_type_from_symmetry<LatticeSetting::primitive>(
    SymmetryOperations const &, Matrix3d const &, double);

} // namespace cppcrystal::spacegroup
