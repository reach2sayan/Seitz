#include <spglib/symmetry/pointgroup.hpp>

#include <array>
#include <cstdlib>
#include <optional>

// Port of pointgroup.c (3D space-group path). Two static tables drive it:
//   * kPointgroupData: for each of the 32 point groups, the histogram of
//     rotation types {-6,-4,-3,-2,-1,1,2,3,4,6} plus its
//     symbols/holohedry/Laue.
//   * kRotAxes: 73 candidate rotation axes used to pick conventional axes.
// Identification matches the input's rotation-type histogram against the table;
// the transformation matrix selects conventional axes per Laue class.
namespace spglib::symmetry {

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
enum class RotationType {
  rotoinversion_6, // -6
  rotoinversion_4, // -4
  rotoinversion_3, // -3
  mirror,          // -2
  inversion,       // -1
  identity,        //  1
  rotation_2,      //  2
  rotation_3,      //  3
  rotation_4,      //  4
  rotation_6,      //  6
};
constexpr std::size_t kRotationTypeCount = 10;

// Classify a rotation by its determinant and trace; std::nullopt if it is not a
// crystallographic rotation (a defensive case that valid point groups avoid).
[[nodiscard]] std::optional<RotationType>
rotation_type(Matrix3i const &rot) noexcept {
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

// De-duplicate rotations by value (ptg_get_pointsymmetry).
[[nodiscard]] PointSymmetry
unique_rotations(std::span<Matrix3i const> rotations) {
  PointSymmetry out;
  for (Matrix3i const &r : rotations) {
    bool seen = false;
    for (Matrix3i const &u : out)
      if (u == r) {
        seen = true;
        break;
      }
    if (!seen && out.size() < out.capacity())
      out.push_back(r);
  }
  return out;
}

[[nodiscard]] bool class_table(std::array<int, kRotationTypeCount> &table,
                               PointSymmetry const &ps) {
  table.fill(0);
  for (Matrix3i const &r : ps) {
    auto const t = rotation_type(r);
    if (!t)
      return false;
    ++table[static_cast<std::size_t>(*t)];
  }
  return true;
}

[[nodiscard]] int pointgroup_number(PointSymmetry const &ps) {
  std::array<int, kRotationTypeCount> table{};
  if (!class_table(table, ps))
    return 0;
  for (int i = 1; i < 33; ++i)
    if (kPointgroupData[static_cast<std::size_t>(i)].table == table)
      return i;
  return 0;
}

[[nodiscard]] Matrix3i proper_rotation(Matrix3i const &rot) {
  return rot.determinant() == -1 ? Matrix3i(-rot) : rot;
}

// Index into kRotAxes that is the rotation axis (eigenvector) of a proper
// rotation; std::nullopt for the identity / when none matches
// (get_rotation_axis).
[[nodiscard]] std::optional<int> rotation_axis(Matrix3i const &proper_rot) {
  if (proper_rot == kIdentity)
    return std::nullopt;
  for (int i = 0; i < kNumRotAxes; ++i) {
    Vector3i const v = proper_rot * rot_axis(i);
    if (v == rot_axis(i))
      return i;
  }
  return std::nullopt;
}

// Indices of kRotAxes orthogonal to the axis of proper_rot
// (get_orthogonal_axis).
[[nodiscard]] std::vector<int> orthogonal_axes(Matrix3i const &proper_rot,
                                               int rot_order) {
  Matrix3i sum_rot = kIdentity;
  Matrix3i rot = kIdentity;
  for (int i = 0; i < rot_order - 1; ++i) {
    rot = proper_rot * rot;
    sum_rot += rot;
  }
  std::vector<int> result;
  for (int i = 0; i < kNumRotAxes; ++i)
    if ((sum_rot * rot_axis(i)).isZero())
      result.push_back(i);
  return result;
}

// 1 if axis_vec == kRotAxes[idx], -1 if == -kRotAxes[idx], else 0.
[[nodiscard]] int axis_sign(Vector3i const &axis_vec, int idx) {
  if (axis_vec == rot_axis(idx))
    return 1;
  if (axis_vec == Vector3i(-rot_axis(idx)))
    return -1;
  return 0;
}

// A chosen conventional axis: an index into kRotAxes plus a direction sign.
// Replaces spglib's convention of encoding a negated axis as index +
// kNumRotAxes.
struct SignedAxis {
  int index = 0; // 0..kNumRotAxes-1
  int sign = 1;  // +1 or -1
};
using AxisTriple = std::array<SignedAxis, 3>;

// Build the integer transformation matrix whose columns are the chosen axes
// (set_transformation_matrix).
[[nodiscard]] Matrix3i transformation_from_axes(AxisTriple const &axes) {
  Matrix3i tmat;
  for (int j = 0; j < 3; ++j) {
    auto const &a = axes[static_cast<std::size_t>(j)];
    tmat.col(j) = a.sign * rot_axis(a.index);
  }
  return tmat;
}

bool laue2m(AxisTriple &axes, PointSymmetry const &ps) {
  Matrix3i prop_rot = kIdentity;
  std::optional<int> two_fold;
  for (Matrix3i const &r : ps) {
    prop_rot = proper_rotation(r);
    if (prop_rot.trace() == -1) { // two-fold
      two_fold = rotation_axis(prop_rot);
      break;
    }
  }
  if (!two_fold)
    return false;
  axes[1] = {*two_fold, 1};

  auto const ortho = orthogonal_axes(prop_rot, 2);
  if (ortho.empty())
    return false;

  int min_norm = 8;
  bool found = false;
  for (int idx : ortho) {
    int const norm = rot_axis(idx).squaredNorm();
    if (norm < min_norm) {
      min_norm = norm;
      axes[0] = {idx, 1};
      found = true;
    }
  }
  if (!found)
    return false;

  min_norm = 8;
  found = false;
  for (int idx : ortho) {
    int const norm = rot_axis(idx).squaredNorm();
    if (norm < min_norm && idx != axes[0].index) {
      min_norm = norm;
      axes[2] = {idx, 1};
      found = true;
    }
  }
  return found;
}

// Laue classes 4/m, 4/mmm, -3, -3m, 6/m, 6/mmm (laue_one_axis).
bool laue_one_axis(AxisTriple &axes, PointSymmetry const &ps, int rot_order) {
  Matrix3i prop_rot = kIdentity;
  std::optional<int> principal;
  for (Matrix3i const &r : ps) {
    prop_rot = proper_rotation(r);
    if ((rot_order == 4 && prop_rot.trace() == 1) ||
        (rot_order == 3 && prop_rot.trace() == 0)) {
      principal = rotation_axis(prop_rot);
      break;
    }
  }
  if (!principal)
    return false;
  axes[2] = {*principal, 1};

  auto const ortho = orthogonal_axes(prop_rot, rot_order);
  if (ortho.empty())
    return false;

  AxisTriple tmp;
  tmp[2] = axes[2];
  bool resolved = false;
  for (int first : ortho) {
    tmp[0] = {first, 1};
    Vector3i const axis_vec = prop_rot * rot_axis(first);
    int second = 0;
    int sign = 0;
    for (int cand : ortho) {
      sign = axis_sign(axis_vec, cand);
      if (sign != 0) {
        second = cand;
        break;
      }
    }
    if (sign == 0)
      continue; // no matching second axis
    tmp[1] = {second, sign};
    if (std::abs(transformation_from_axes(tmp).determinant()) < 4) {
      // < 4 avoids the F-centred (det 4) choice.
      axes[0] = tmp[0];
      axes[1] = tmp[1];
      resolved = true;
      break;
    }
  }
  if (!resolved)
    return false;

  if (transformation_from_axes(axes).determinant() < 0)
    std::swap(axes[0], axes[1]);
  return true;
}

void sort_axes(AxisTriple &axes) {
  if (axes[1].index > axes[2].index)
    std::swap(axes[1], axes[2]);
  if (axes[0].index > axes[1].index)
    std::swap(axes[0], axes[1]);
  if (axes[1].index > axes[2].index)
    std::swap(axes[1], axes[2]);
  if (transformation_from_axes(axes).determinant() < 0)
    std::swap(axes[1], axes[2]);
}

// Laue classes mmm, m-3, m-3m (lauennn, 3D path).
bool lauennn(AxisTriple &axes, PointSymmetry const &ps, int rot_order) {
  std::array<int, 3> idx{};
  int count = 0;
  for (Matrix3i const &r : ps) {
    Matrix3i const prop_rot = proper_rotation(r);
    if ((prop_rot.trace() == -1 && rot_order == 2) ||
        (prop_rot.trace() == 1 && rot_order == 4)) {
      auto const axis = rotation_axis(prop_rot);
      if (!axis)
        continue;
      bool dup = false;
      for (int i = 0; i < count; ++i)
        if (idx[static_cast<std::size_t>(i)] == *axis)
          dup = true;
      if (!dup && count < 3)
        idx[static_cast<std::size_t>(count++)] = *axis;
    }
  }
  if (count != 3)
    return false; // mmm / m-3 / m-3m each have exactly three such axes
  for (int i = 0; i < 3; ++i)
    axes[static_cast<std::size_t>(i)] = {idx[static_cast<std::size_t>(i)], 1};
  sort_axes(axes);
  return true;
}

[[nodiscard]] bool get_axes(AxisTriple &axes, Laue laue,
                            PointSymmetry const &ps) {
  switch (laue) {
  case Laue::laue_1:
    axes = {SignedAxis{0, 1}, SignedAxis{1, 1}, SignedAxis{2, 1}};
    return true;
  case Laue::laue_2m:
    return laue2m(axes, ps);
  case Laue::laue_mmm:
    return lauennn(axes, ps, 2);
  case Laue::laue_4m:
  case Laue::laue_4mmm:
    return laue_one_axis(axes, ps, 4);
  case Laue::laue_3:
  case Laue::laue_3m:
  case Laue::laue_6m:
  case Laue::laue_6mmm:
    return laue_one_axis(axes, ps, 3);
  case Laue::laue_m3:
    return lauennn(axes, ps, 2);
  case Laue::laue_m3m:
    return lauennn(axes, ps, 4);
  default:
    return false;
  }
}

} // namespace

PointGroup pointgroup_by_number(int number) noexcept {
  if (number < 1 || number > 32)
    return {};
  PgEntry const &e = kPointgroupData[static_cast<std::size_t>(number)];
  return {number, e.symbol, e.schoenflies, e.holohedry, e.laue};
}

int identify_pointgroup_number(std::span<Matrix3i const> rotations) noexcept {
  return pointgroup_number(unique_rotations(rotations));
}

Result<PointgroupTransform>
get_pointgroup(std::span<Matrix3i const> rotations) {
  PointSymmetry const ps = unique_rotations(rotations);
  int const pg_num = pointgroup_number(ps);
  if (pg_num == 0)
    return leaf::new_error(e_pointgroup_not_found{});

  PointgroupTransform result;
  result.pointgroup = pointgroup_by_number(pg_num);
  AxisTriple axes;
  if (!get_axes(axes, result.pointgroup.laue, ps))
    return leaf::new_error(e_pointgroup_not_found{});
  result.transformation = transformation_from_axes(axes);
  return result;
}

std::vector<Matrix3i> rotations_of(SymmetryOperations const &ops) {
  std::vector<Matrix3i> r;
  r.reserve(ops.size());
  for (auto const &op : ops)
    r.push_back(op.rotation);
  return r;
}

} // namespace spglib::symmetry
