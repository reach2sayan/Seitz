#include <cppcrystal/symmetry/pointgroup.hpp>

#include <cppcrystal/core/matrix_order.hpp>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <optional>
#include <ranges>
#include <utility>

// Two static tables drive point-group identification:
//   * kPointgroupData: for each of the 32 point groups, the histogram of
//     rotation types {-6,-4,-3,-2,-1,1,2,3,4,6} plus its
//     symbols/holohedry/Laue.
//   * kRotAxes: 73 candidate rotation axes used to pick conventional axes.
// Identification matches the input's rotation-type histogram against the table;
// the transformation matrix selects conventional axes per Laue class.
namespace cppcrystal::symmetry {

namespace {

struct PgEntry {
  std::array<int, 10> table;
  std::string_view symbol;
  std::string_view schoenflies;
  Holohedry holohedry;
  Laue laue;
};

constexpr std::array<PgEntry, 33> kPointgroupData = {{
    {{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "", "", Holohedry::none, Laue::none},
    {{0, 0, 0, 0, 0, 1, 0, 0, 0, 0},
     "1",
     "C1",
     Holohedry::triclinic,
     Laue::laue_1},
    {{0, 0, 0, 0, 1, 1, 0, 0, 0, 0},
     "-1",
     "Ci",
     Holohedry::triclinic,
     Laue::laue_1},
    {{0, 0, 0, 0, 0, 1, 1, 0, 0, 0},
     "2",
     "C2",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 1, 0, 1, 0, 0, 0, 0},
     "m",
     "Cs",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
     "2/m",
     "C2h",
     Holohedry::monoclinic,
     Laue::laue_2m},
    {{0, 0, 0, 0, 0, 1, 3, 0, 0, 0},
     "222",
     "D2",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 2, 0, 1, 1, 0, 0, 0},
     "mm2",
     "C2v",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 3, 1, 1, 3, 0, 0, 0},
     "mmm",
     "D2h",
     Holohedry::orthorhombic,
     Laue::laue_mmm},
    {{0, 0, 0, 0, 0, 1, 1, 0, 2, 0},
     "4",
     "C4",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 2, 0, 0, 0, 1, 1, 0, 0, 0},
     "-4",
     "S4",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 2, 0, 1, 1, 1, 1, 0, 2, 0},
     "4/m",
     "C4h",
     Holohedry::tetragonal,
     Laue::laue_4m},
    {{0, 0, 0, 0, 0, 1, 5, 0, 2, 0},
     "422",
     "D4",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 0, 0, 4, 0, 1, 1, 0, 2, 0},
     "4mm",
     "C4v",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 2, 0, 2, 0, 1, 3, 0, 0, 0},
     "-42m",
     "D2d",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 2, 0, 5, 1, 1, 5, 0, 2, 0},
     "4/mmm",
     "D4h",
     Holohedry::tetragonal,
     Laue::laue_4mmm},
    {{0, 0, 0, 0, 0, 1, 0, 2, 0, 0},
     "3",
     "C3",
     Holohedry::trigonal,
     Laue::laue_3},
    {{0, 0, 2, 0, 1, 1, 0, 2, 0, 0},
     "-3",
     "C3i",
     Holohedry::trigonal,
     Laue::laue_3},
    {{0, 0, 0, 0, 0, 1, 3, 2, 0, 0},
     "32",
     "D3",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 0, 3, 0, 1, 0, 2, 0, 0},
     "3m",
     "C3v",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 2, 3, 1, 1, 3, 2, 0, 0},
     "-3m",
     "D3d",
     Holohedry::trigonal,
     Laue::laue_3m},
    {{0, 0, 0, 0, 0, 1, 1, 2, 0, 2},
     "6",
     "C6",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{2, 0, 0, 1, 0, 1, 0, 2, 0, 0},
     "-6",
     "C3h",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{2, 0, 2, 1, 1, 1, 1, 2, 0, 2},
     "6/m",
     "C6h",
     Holohedry::hexagonal,
     Laue::laue_6m},
    {{0, 0, 0, 0, 0, 1, 7, 2, 0, 2},
     "622",
     "D6",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{0, 0, 0, 6, 0, 1, 1, 2, 0, 2},
     "6mm",
     "C6v",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{2, 0, 0, 4, 0, 1, 3, 2, 0, 0},
     "-6m2",
     "D3h",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{2, 0, 2, 7, 1, 1, 7, 2, 0, 2},
     "6/mmm",
     "D6h",
     Holohedry::hexagonal,
     Laue::laue_6mmm},
    {{0, 0, 0, 0, 0, 1, 3, 8, 0, 0},
     "23",
     "T",
     Holohedry::cubic,
     Laue::laue_m3},
    {{0, 0, 8, 3, 1, 1, 3, 8, 0, 0},
     "m-3",
     "Th",
     Holohedry::cubic,
     Laue::laue_m3},
    {{0, 0, 0, 0, 0, 1, 9, 8, 6, 0},
     "432",
     "O",
     Holohedry::cubic,
     Laue::laue_m3m},
    {{0, 6, 0, 6, 0, 1, 3, 8, 0, 0},
     "-43m",
     "Td",
     Holohedry::cubic,
     Laue::laue_m3m},
    {{0, 6, 8, 9, 1, 1, 9, 8, 6, 0},
     "m-3m",
     "Oh",
     Holohedry::cubic,
     Laue::laue_m3m},
}};

constexpr int kNumRotAxes = 73;
constexpr int kRotAxes[kNumRotAxes][3] = {
    {1, 0, 0},   {0, 1, 0},   {0, 0, 1},   {0, 1, 1},   {1, 0, 1},
    {1, 1, 0},   {0, 1, -1},  {-1, 0, 1},  {1, -1, 0},  {1, 1, 1},
    {-1, 1, 1},  {1, -1, 1},  {1, 1, -1},  {0, 1, 2},   {2, 0, 1},
    {1, 2, 0},   {0, 2, 1},   {1, 0, 2},   {2, 1, 0},   {0, -1, 2},
    {2, 0, -1},  {-1, 2, 0},  {0, -2, 1},  {1, 0, -2},  {-2, 1, 0},
    {2, 1, 1},   {1, 2, 1},   {1, 1, 2},   {2, -1, -1}, {-1, 2, -1},
    {-1, -1, 2}, {2, 1, -1},  {-1, 2, 1},  {1, -1, 2},  {2, -1, 1},
    {1, 2, -1},  {-1, 1, 2},  {3, 1, 2},   {2, 3, 1},   {1, 2, 3},
    {3, 2, 1},   {1, 3, 2},   {2, 1, 3},   {3, -1, 2},  {2, 3, -1},
    {-1, 2, 3},  {3, -2, 1},  {1, 3, -2},  {-2, 1, 3},  {3, -1, -2},
    {-2, 3, -1}, {-1, -2, 3}, {3, -2, -1}, {-1, 3, -2}, {-2, -1, 3},
    {3, 1, -2},  {-2, 3, 1},  {1, -2, 3},  {3, 2, -1},  {-1, 3, 2},
    {2, -1, 3},  {1, 1, 3},   {-1, 1, 3},  {1, -1, 3},  {-1, -1, 3},
    {1, 3, 1},   {-1, 3, 1},  {1, 3, -1},  {-1, 3, -1}, {3, 1, 1},
    {3, 1, -1},  {3, -1, 1},  {3, -1, -1},
};

[[nodiscard]] Vector3i rot_axis(int i) {
  return {kRotAxes[i][0], kRotAxes[i][1], kRotAxes[i][2]};
}

Matrix3i const kIdentity = Matrix3i::Identity();

// The ten crystallographic rotation types, ordered as the point-group
// class-count table indexes them: {-6,-4,-3,-2,-1,1,2,3,4,6}.
constexpr std::size_t kRotationTypeCount = 10;

// rotation_type() is now public (declared in the header) and defined below,
// after this anonymous namespace; class_table() calls it through that
// declaration.

// De-duplicate rotations by value; distinct rotations beyond the capacity of
// PointSymmetry are dropped, as before.
[[nodiscard]] PointSymmetry
unique_rotations(std::span<Matrix3i const> rotations) {
  PointSymmetry out;
  for (const Matrix3i &r : unique_by_rotation(rotations)) {
    if (out.size() < out.capacity()) {
      out.push_back(r);
    }
  }
  return out;
}

[[nodiscard]] bool class_table(std::array<int, kRotationTypeCount> &table,
                               PointSymmetry const &ps) {
  table.fill(0);
  for (const auto &r : ps) {
    if (const auto t = rotation_type(r)) {
      ++table[static_cast<std::size_t>(*t)];
    } else {
      return false;
    }
  }
  return true;
}

[[nodiscard]] int pointgroup_number(PointSymmetry const &ps) {
  std::array<int, kRotationTypeCount> table{};
  if (!class_table(table, ps)) {
    return 0;
  }
  auto const it = std::ranges::find(kPointgroupData, table, &PgEntry::table);
  return it == kPointgroupData.end()
             ? 0
             : static_cast<int>(it - kPointgroupData.begin());
}

[[nodiscard]] inline Matrix3i proper_rotation(Matrix3i const &rot) {
  return rot.determinant() == -1 ? Matrix3i(-rot) : rot;
}

// Index into kRotAxes that is the rotation axis (eigenvector) of a proper
// rotation; std::nullopt for the identity / when none matches.
[[nodiscard]] std::optional<int> axis_index(Matrix3i const &proper_rot) {
  if (proper_rot != kIdentity) {
    for (int i = 0; i < kNumRotAxes; ++i) {
      Vector3i const v = proper_rot * rot_axis(i);
      if (v == rot_axis(i)) {
        return i;
      }
    }
  }
  return std::nullopt;
}

// Indices of kRotAxes orthogonal to the axis of proper_rot.
[[nodiscard]] std::vector<int> orthogonal_axes(Matrix3i const &proper_rot,
                                               int rot_order) {
  Matrix3i sum_rot = kIdentity;
  Matrix3i rot = kIdentity;
  for (int i = 0; i < rot_order - 1; ++i) {
    rot = proper_rot * rot;
    sum_rot += rot;
  }
  std::vector<int> result;
  result.reserve(kNumRotAxes);
  std::ranges::copy_if(std::views::iota(0, kNumRotAxes),
                       std::back_inserter(result),
                       [&](int i) { return (sum_rot * rot_axis(i)).isZero(); });
  return result;
}

// 1 if axis_vec == kRotAxes[idx], -1 if == -kRotAxes[idx], else 0.
[[nodiscard]] int axis_sign(Vector3i const &axis_vec, int idx) {
  if (axis_vec == rot_axis(idx)) {
    return 1;
  } else if (axis_vec == Vector3i(-rot_axis(idx))) {
    return -1;
  }
  return 0;
}

// A chosen conventional axis: an index into kRotAxes plus a direction sign
// (rather than encoding a negated axis as index + kNumRotAxes).
struct SignedAxis {
  int index = 0; // 0..kNumRotAxes-1
  int sign = 1;  // +1 or -1
};
using AxisTriple = std::array<SignedAxis, 3>;

// Build the integer transformation matrix whose columns are the chosen axes.
[[nodiscard]] Matrix3i transformation_from_axes(AxisTriple const &axes) {
  Matrix3i tmat;
  for (auto const &[j, a] : axes | std::views::enumerate) {
    tmat.col(j) = a.sign * rot_axis(a.index);
  }
  return tmat;
}

// Squared length of candidate axis kRotAxes[idx].
[[nodiscard]] int axis_sqnorm(int idx) { return rot_axis(idx).squaredNorm(); }

// The shortest axis (index into kRotAxes) of `range`; min_element keeps the
// first occurrence on ties (the index-order tie-break). nullopt when empty.
template <std::ranges::forward_range R>
[[nodiscard]] std::optional<int> shortest_axis(R &&range) {
  auto const it = std::ranges::min_element(range, {}, axis_sqnorm);
  if (it == std::ranges::end(range)) {
    return std::nullopt;
  }
  return *it;
}

// The distinguished axis of a Laue class: the first proper rotation of
// rot_order (2, 3 or 4, identified by its trace) paired with its rotation-axis
// index. nullopt if no such operation (or no axis for it) exists.
struct Principal {
  Matrix3i prop_rot;
  int axis;
};
[[nodiscard]] std::optional<Principal> principal_axis(PointSymmetry const &ps,
                                                      int rot_order) {
  for (Matrix3i const &r : ps) {
    Matrix3i const prop_rot = proper_rotation(r);
    if ((rot_order == 4 && prop_rot.trace() == 1) ||
        (rot_order == 3 && prop_rot.trace() == 0) ||
        (rot_order == 2 && prop_rot.trace() == -1)) {
      return axis_index(prop_rot).transform(
          [&](int axis) { return Principal{prop_rot, axis}; });
    }
  }
  return std::nullopt;
}

bool laue2m(AxisTriple &axes, PointSymmetry const &ps) {
  auto const two_fold = principal_axis(ps, 2);
  if (!two_fold) {
    return false;
  }
  axes[1] = {two_fold->axis, 1};

  auto const ortho = orthogonal_axes(two_fold->prop_rot, 2);
  if (ortho.empty()) {
    return false;
  }

  // a and b: the two shortest orthogonal axes (distinct).
  axes[0] = {*shortest_axis(ortho), 1};

  auto const second = shortest_axis(ortho | std::views::filter([&](int idx) {
                                      return idx != axes[0].index;
                                    }));
  if (!second) {
    return false;
  }
  axes[2] = {*second, 1};
  return true;
}

// The two in-plane axes (a, b) for a one-axis Laue class: the first orthogonal
// axis whose image under prop_rot is, up to sign, another orthogonal axis,
// chosen so the cell is not F-centred (|det| < 4). nullopt if none qualifies
// (also covers an empty orthogonal-axis set).
[[nodiscard]] std::optional<std::pair<SignedAxis, SignedAxis>>
in_plane_axes(Principal const &p, std::vector<int> const &ortho) {
  for (int const first : ortho) {
    Vector3i const axis_vec = p.prop_rot * rot_axis(first);
    for (int const cand : ortho) {
      int const sign = axis_sign(axis_vec, cand);
      if (sign == 0) {
        continue;
      }
      AxisTriple const tmp{{SignedAxis{first, 1}, SignedAxis{cand, sign},
                            SignedAxis{p.axis, 1}}};
      if (std::abs(transformation_from_axes(tmp).determinant()) < 4) {
        return std::pair{SignedAxis{first, 1}, SignedAxis{cand, sign}};
      }
      break; // only the first matching second axis is considered per `first`
    }
  }
  return std::nullopt;
}

// Laue classes 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm.
[[nodiscard]] std::optional<AxisTriple> laue_one_axis(PointSymmetry const &ps,
                                                      int rot_order) {
  return principal_axis(ps, rot_order).and_then([&](Principal const &p) {
    return in_plane_axes(p, orthogonal_axes(p.prop_rot, rot_order))
        .transform([&](std::pair<SignedAxis, SignedAxis> const &ab) {
          AxisTriple axes{{ab.first, ab.second, SignedAxis{p.axis, 1}}};
          if (transformation_from_axes(axes).determinant() < 0) {
            std::swap(axes[0], axes[1]);
          }
          return axes;
        });
  });
}

void sort_axes(AxisTriple &axes) {
  std::ranges::sort(axes, {}, &SignedAxis::index);
  if (transformation_from_axes(axes).determinant() < 0) {
    std::swap(axes[1], axes[2]);
  }
}

// The aperiodic-axis component of candidate axis `idx` (rot_axes[idx][ap]).
[[nodiscard]] constexpr int aperiodic_component(int idx, int aperiodic_axis) {
  return kRotAxes[idx][aperiodic_axis];
}

// Of the three chosen axes, exactly two must lie in the periodic plane
// (aperiodic component 0) and one along the aperiodic axis (component +/-1).
// Move the aperiodic axis to c, then orient for a positive determinant. False
// for an invalid (e.g. inclined) configuration.
[[nodiscard]] bool layer_sort_axes(AxisTriple &axes, int aperiodic_axis) {
  auto const component = [&](int i) {
    return aperiodic_component(axes[static_cast<std::size_t>(i)].index,
                               aperiodic_axis);
  };
  auto const in_plane = [&](int i) { return component(i) == 0; };
  auto const along_aperiodic = [&](int i) {
    int const c = component(i);
    return c == 1 || c == -1;
  };
  auto const positions = std::views::iota(0, 3);
  if (std::ranges::count_if(positions, in_plane) != 2 ||
      std::ranges::count_if(positions, along_aperiodic) != 1) {
    return false;
  }
  int const axis_pos = *std::ranges::find_if(positions, along_aperiodic);
  std::swap(axes[static_cast<std::size_t>(axis_pos)], axes[2]);
  if (transformation_from_axes(axes).determinant() < 0) {
    std::swap(axes[0], axes[1]);
  }
  return true;
}

// Layer LAUE2M: unlike the 3D case the two-fold axis becomes axis a (its
// position relative to the aperiodic axis distinguishes oblique vs. rectangular
// monoclinic layers); the remaining two axes are chosen accordingly.
[[nodiscard]] bool layer_laue2m(AxisTriple &axes, PointSymmetry const &ps,
                                int aperiodic_axis) {
  auto const two_fold = principal_axis(ps, 2);
  if (!two_fold) {
    return false;
  }
  axes[0] = {two_fold->axis, 1};

  auto const ortho = orthogonal_axes(two_fold->prop_rot, 2);
  if (ortho.empty()) {
    return false;
  }

  int const a0 = aperiodic_component(two_fold->axis, aperiodic_axis);
  if (a0 == 1 || a0 == -1) {
    // Monoclinic/oblique: the two-fold is along the aperiodic axis; a and b are
    // the two shortest orthogonal axes.
    auto const first = shortest_axis(ortho);
    if (!first) {
      return false;
    }
    auto const second =
        shortest_axis(ortho | std::views::filter(
                                  [&](int idx) { return idx != axes[1].index; }));
    if (!second) {
      return false;
    }
    axes[1] = {*first, 1};
    axes[2] = {*second, 1};
  } else if (a0 == 0) {
    // Monoclinic/rectangular: the second axis lies in the periodic plane, the
    // third along the aperiodic axis.
    auto const in_plane =
        shortest_axis(ortho | std::views::filter([&](int idx) {
                        return aperiodic_component(idx, aperiodic_axis) == 0;
                      }));
    if (!in_plane) {
      return false;
    }

    auto const out_plane =
        shortest_axis(ortho | std::views::filter([&](int idx) {
                        int const c = aperiodic_component(idx, aperiodic_axis);
                        return c == 1 || c == -1;
                      }));
    if (!out_plane) {
      return false;
    }

    axes[1] = {*in_plane, 1};
    axes[2] = {*out_plane, 1};
  } else {
    return false;
  }
  return true;
}

// Laue classes mmm, m-3, m-3m. For a layer cell the three axes are reordered so
// the aperiodic axis is c.
bool lauennn(AxisTriple &axes, PointSymmetry const &ps, int rot_order,
             std::optional<int> aperiodic_axis) {
  std::array<int, 3> idx{};
  int count = 0;
  for (Matrix3i const &r : ps) {
    Matrix3i const prop_rot = proper_rotation(r);
    if ((prop_rot.trace() == -1 && rot_order == 2) ||
        (prop_rot.trace() == 1 && rot_order == 4)) {
      auto const axis = axis_index(prop_rot);
      if (!axis) {
        continue;
      }
      auto const seen = std::ranges::subrange(idx.begin(), idx.begin() + count);
      if (std::ranges::find(seen, *axis) == seen.end() && count < 3) {
        idx[static_cast<std::size_t>(count++)] = *axis;
      }
    }
  }
  if (count != 3) {
    return false; // mmm / m-3 / m-3m each have exactly three such axes
  }
  for (auto [a, id] : std::views::zip(axes, idx)) {
    a = {id, 1};
  }
  if (aperiodic_axis) {
    return layer_sort_axes(axes, *aperiodic_axis);
  }
  sort_axes(axes);
  return true;
}

[[nodiscard]] bool get_axes(AxisTriple &axes, Laue laue,
                            PointSymmetry const &ps,
                            std::optional<int> aperiodic_axis) {
  auto const assign_axes = [&](int rot_order) {
    return laue_one_axis(ps, rot_order)
        .transform([&](AxisTriple const &a) {
          axes = a;
          return true;
        })
        .value_or(false);
  };
  switch (laue) {
  case Laue::laue_1:
    axes = {SignedAxis{0, 1}, SignedAxis{1, 1}, SignedAxis{2, 1}};
    return true;
  case Laue::laue_2m:
    return aperiodic_axis ? layer_laue2m(axes, ps, *aperiodic_axis)
                          : laue2m(axes, ps);
  case Laue::laue_mmm:
    return lauennn(axes, ps, 2, aperiodic_axis);
  case Laue::laue_4m:
  case Laue::laue_4mmm:
    // The 4-fold axis is the aperiodic axis for a layer; laue_one_axis already
    // places it at c, so no aperiodic-specific handling is needed.
    return assign_axes(4);
  case Laue::laue_3:
  case Laue::laue_3m:
  case Laue::laue_6m:
  case Laue::laue_6mmm:
    return assign_axes(3);
  case Laue::laue_m3:
    return lauennn(axes, ps, 2, aperiodic_axis);
  case Laue::laue_m3m:
    return lauennn(axes, ps, 4, aperiodic_axis);
  default:
    return false;
  }
}

} // namespace

std::optional<RotationType> rotation_type(Matrix3i const &rot) noexcept {
  int const det = rot.determinant();
  int const tr = rot.trace();
  if (det == -1) {
    switch (tr) {
    case -2:
      return RotationType::rotoinversion_6;
    case -1:
      return RotationType::rotoinversion_4;
    case 0:
      return RotationType::rotoinversion_3;
    case 1:
      return RotationType::mirror;
    case -3:
      return RotationType::inversion;
    default:
      return std::nullopt;
    }
  }
  switch (tr) {
  case 3:
    return RotationType::identity;
  case -1:
    return RotationType::rotation_2;
  case 0:
    return RotationType::rotation_3;
  case 1:
    return RotationType::rotation_4;
  case 2:
    return RotationType::rotation_6;
  default:
    return std::nullopt;
  }
}

int rotation_order(Matrix3i const &rot) noexcept {
  // Signed orders in RotationType declaration order.
  constexpr std::array<int, 10> kOrders{-6, -4, -3, -2, -1, 1, 2, 3, 4, 6};
  return rotation_type(rot)
      .transform([&](RotationType t) {
        return kOrders[static_cast<std::size_t>(std::to_underlying(t))];
      })
      .value_or(0);
}

std::optional<Vector3i> rotation_axis(Matrix3i const &rot) {
  auto const idx = axis_index(proper_rotation(rot));
  return idx ? std::optional<Vector3i>{rot_axis(*idx)} : std::nullopt;
}

PointGroup pointgroup_by_number(int number) noexcept {
  if (number < 1 || number > 32) {
    return {};
  }
  PgEntry const &e = kPointgroupData[static_cast<std::size_t>(number)];
  // CrystalClass enumerators are aligned to the point-group numbering (1..32).
  return {number, e.symbol, e.schoenflies, e.holohedry, e.laue,
          static_cast<CrystalClass>(number)};
}

int identify_pointgroup_number(std::span<Matrix3i const> rotations) noexcept {
  return pointgroup_number(unique_rotations(rotations));
}

Result<PointgroupTransform> get_pointgroup(std::span<Matrix3i const> rotations,
                                           std::optional<int> aperiodic_axis) {
  PointSymmetry const ps = unique_rotations(rotations);
  int const pg_num = pointgroup_number(ps);
  if (pg_num == 0) {
    return leaf::new_error(e_pointgroup_not_found{});
  }
  // Layer groups have no cubic point groups (numbers 28..32).
  if (aperiodic_axis && pg_num >= 28) {
    return leaf::new_error(e_pointgroup_not_found{});
  }

  PointgroupTransform result;
  result.pointgroup = pointgroup_by_number(pg_num);
  AxisTriple axes;
  if (!get_axes(axes, result.pointgroup.laue, ps, aperiodic_axis)) {
    return leaf::new_error(e_pointgroup_not_found{});
  }
  result.transformation = transformation_from_axes(axes);
  return result;
}

std::vector<Matrix3i> rotations_of(SymmetryOperations const &ops) {
  std::vector<Matrix3i> r;
  r.reserve(ops.size());

  std::ranges::transform(ops, std::back_inserter(r),
                         &SymmetryOperation::rotation);

  return r;
}

} // namespace cppcrystal::symmetry
