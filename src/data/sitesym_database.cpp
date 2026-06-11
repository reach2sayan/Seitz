#include <spglib/data/sitesym_database.hpp>

#include <spglib/data/sitesym_tables.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <ranges>
#include <string_view>

namespace spglib::data {

namespace {

constexpr int kRowBase = 45;
constexpr int kRotEncMod = kRowBase * kRowBase * kRowBase;
constexpr int kRow2 = kRowBase * kRowBase;

// Compile-time decoded form of one Wyckoff representative: rotation as 9
// row-major entries in {-2..2}/{-1..1}, translation as 3 base-24 numerators
// (value n/24), plus the multiplicity. The whole encoded table is unpacked once
// at compile time so wyckoff_coordinate only assembles Eigen at runtime, never
// decodes. Mirrors the kDecodedOps pattern in spg_database.cpp.
struct WyckoffLit {
  std::array<std::int8_t, 9> rot;
  std::array<std::int8_t, 3> trans_num;
  int multiplicity;
};

constexpr auto kWyckoffDecoded = [] {
  std::array<WyckoffLit, kCoordinatesFirst.size()> out{};
  for (std::size_t k = 0; k < kCoordinatesFirst.size(); ++k) {
    int const code = kCoordinatesFirst[k];

    // Rotation: base-45 row codes, each digit base-9/3 decoded to
    // {-2..2}/{-1..1}.
    int const rot_enc = code % kRotEncMod;
    std::array<int, 3> const rows = {
        rot_enc / kRow2, (rot_enc % kRow2) / kRowBase, rot_enc % kRowBase};
    for (std::size_t i = 0; i < 3; ++i) {
      out[k].rot[i * 3 + 0] = static_cast<std::int8_t>(rows[i] / 9 - 2);
      out[k].rot[i * 3 + 1] = static_cast<std::int8_t>((rows[i] % 9) / 3 - 1);
      out[k].rot[i * 3 + 2] = static_cast<std::int8_t>(rows[i] % 3 - 1);
    }

    // Translation: base-24 per axis (integer digits), each a fraction of 24.
    int const trans_enc = code / kRotEncMod;
    out[k].trans_num[0] = static_cast<std::int8_t>(trans_enc / 576);
    out[k].trans_num[1] = static_cast<std::int8_t>((trans_enc % 576) / 24);
    out[k].trans_num[2] = static_cast<std::int8_t>(trans_enc % 24);

    out[k].multiplicity = kMultiplicities[k];
  }
  return out;
}();

// Index 1 is the general position's identity representative x' = E.x (zero
// translation): a compile-time check that the decode math is intact.
static_assert(kWyckoffDecoded[1].rot ==
              std::array<std::int8_t, 9>{1, 0, 0, 0, 1, 0, 0, 0, 1});
static_assert(kWyckoffDecoded[1].trans_num ==
              std::array<std::int8_t, 3>{0, 0, 0});

// Per-Hall Wyckoff ranges, precomputed once at compile time from
// kPositionWyckoff (the hall_number domain is fixed at 1..530). Index 0 is the
// out-of-range sentinel ({0, 0}); valid Hall numbers index directly.
// kPositionWyckoff has kNumHalls + 2 entries (0..531), so this table holds
// kNumHalls + 1.
constexpr auto kWyckoffRanges = [] {
  std::array<WyckoffRange, kPositionWyckoff.size() - 1> t{};
  for (std::size_t const h : std::views::iota(std::size_t{1}, t.size())) {
    int const start = kPositionWyckoff[h];
    t[h] = WyckoffRange{start, kPositionWyckoff[h + 1] - start};
  }
  return t;
}();
} // namespace

WyckoffCoordinate wyckoff_coordinate(int index) {
  WyckoffLit const &d = kWyckoffDecoded[static_cast<std::size_t>(index)];

  Matrix3i rot;
  rot << d.rot[0], d.rot[1], d.rot[2], d.rot[3], d.rot[4], d.rot[5], d.rot[6],
      d.rot[7], d.rot[8];

  Vector3d const trans{d.trans_num[0] / 24.0, d.trans_num[1] / 24.0,
                       d.trans_num[2] / 24.0};

  return {rot, trans, d.multiplicity};
}

WyckoffRange wyckoff_indices(int hall_number) {
  bool const in_range =
      hall_number >= 1 &&
      static_cast<std::size_t>(hall_number) < kWyckoffRanges.size();
  return kWyckoffRanges[static_cast<std::size_t>(in_range ? hall_number : 0)];
}

std::string_view site_symmetry_symbol(int index) {
  return kSiteSymmetrySymbols[static_cast<std::size_t>(index)];
}

} // namespace spglib::data
