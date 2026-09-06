#include "symmetry/pointgroup.hpp"
#include <seitz/core/operation_set.hpp>

#include "core/matrix_order.hpp"

#include <boost/container/static_vector.hpp>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <numeric>
#include <optional>
#include <ranges>
#include <utility>

// Two static tables drive identification: kPointgroupData holds each of the 32
// groups' rotation-type histogram {-6,-4,-3,-2,-1,1,2,3,4,6} and its
// symbols/holohedry/Laue; kRotAxes holds 73 candidate axes. The input's
// histogram picks the group, then the Laue class picks conventional axes.
namespace seitz::symmetry {

namespace {

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

// Crystallographic type of a single rotation, in a lattice basis, ordered as
// the class-count table indexes them.
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

// A crystallographic rotation is fixed by its determinant (+-1) and its trace
// (-3..3), so the classification is a 2 x 7 lookup rather than a pair of
// switches. Built once at compile time; unset cells are non-crystallographic.
inline constexpr auto kTypeByDetTrace = [] {
  std::array<std::array<std::optional<RotationType>, 7>, 2> t{};
  auto at = [&t](int det, int trace) -> std::optional<RotationType> & {
    return t[det == -1 ? 0 : 1][static_cast<std::size_t>(trace + 3)];
  };
  at(-1, -2) = RotationType::rotoinversion_6;
  at(-1, -1) = RotationType::rotoinversion_4;
  at(-1, 0) = RotationType::rotoinversion_3;
  at(-1, 1) = RotationType::mirror;
  at(-1, -3) = RotationType::inversion;
  at(1, 3) = RotationType::identity;
  at(1, -1) = RotationType::rotation_2;
  at(1, 0) = RotationType::rotation_3;
  at(1, 1) = RotationType::rotation_4;
  at(1, 2) = RotationType::rotation_6;
  return t;
}();

// Classify a single integer rotation by determinant and trace; std::nullopt if
// it is not a crystallographic rotation.
[[nodiscard]] constexpr std::optional<RotationType>
rotation_type(int det, int trace) noexcept {
  if ((det != 1 && det != -1) || trace < -3 || trace > 3) {
    return std::nullopt;
  }
  return kTypeByDetTrace[det == -1 ? 0 : 1]
                        [static_cast<std::size_t>(trace + 3)];
}

[[nodiscard]] std::optional<RotationType>
rotation_type(Matrix3i const &rot) noexcept {
  return rotation_type(rot.determinant(), rot.trace());
}

// The ten crystallographic types are exactly the reachable (det, trace) cells.
static_assert(rotation_type(1, 3) == RotationType::identity);
static_assert(rotation_type(-1, -3) == RotationType::inversion);
static_assert(rotation_type(-1, 1) == RotationType::mirror);
static_assert(!rotation_type(1, -2).has_value()); // no proper 6-bar
static_assert(!rotation_type(0, 0).has_value());  // singular
static_assert(!rotation_type(2, 3).has_value());  // not unimodular

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

// Rotation-type histogram of a point symmetry; nullopt if any rotation is not
// crystallographic.
[[nodiscard]] std::optional<std::array<int, kRotationTypeCount>>
class_table(PointSymmetry const &ps) {
  std::array<int, kRotationTypeCount> table{};
  for (const auto &r : ps) {
    auto const t = rotation_type(r);
    if (!t) {
      return std::nullopt;
    }
    ++table[static_cast<std::size_t>(*t)];
  }
  return table;
}

[[nodiscard]] int pointgroup_number(PointSymmetry const &ps) {
  auto const table = class_table(ps);
  if (!table) {
    return 0;
  }
  auto const it = std::ranges::find(kPointgroupData, *table, &PgEntry::table);
  return it == kPointgroupData.end()
             ? 0
             : static_cast<int>(it - kPointgroupData.begin());
}

[[nodiscard]] inline Matrix3i proper_rotation(Matrix3i const &rot) {
  return rot.determinant() == -1 ? Matrix3i(-rot) : rot;
}

// kRotAxes as a sorted (axis, index) table for lookup by value. Every axis
// line appears once, as a primitive vector of one sign.
using AxisKey = std::array<int, 3>;
constexpr auto kAxisByVector = [] {
  std::array<std::pair<AxisKey, int>, kNumRotAxes> t{};
  for (int i = 0; i < kNumRotAxes; ++i) {
    t[static_cast<std::size_t>(i)] = {
        {kRotAxes[i][0], kRotAxes[i][1], kRotAxes[i][2]}, i};
  }
  std::ranges::sort(t, {}, &std::pair<AxisKey, int>::first);
  return t;
}();

[[nodiscard]] std::optional<int> table_index(Vector3i const &axis) {
  AxisKey const key{axis[0], axis[1], axis[2]};
  auto const it = std::ranges::lower_bound(kAxisByVector, key, {},
                                           &std::pair<AxisKey, int>::first);
  if (it == kAxisByVector.end() || it->first != key) {
    return std::nullopt;
  }
  return it->second;
}

// The rotation axis of a proper rotation as a primitive integer vector: a
// null vector of R - I (the cross product of two independent rows), reduced
// by its gcd. nullopt for the identity, whose fixed space is everything.
[[nodiscard]] std::optional<Vector3i>
primitive_axis(Matrix3i const &proper_rot) {
  Matrix3i const m = proper_rot - kIdentity;
  for (auto const &[a, b] :
       {std::pair{0, 1}, std::pair{0, 2}, std::pair{1, 2}}) {
    Vector3i const v =
        Vector3i(m.row(a).transpose()).cross(Vector3i(m.row(b).transpose()));
    if (!v.isZero()) {
      int const g =
          std::gcd(std::gcd(std::abs(v[0]), std::abs(v[1])), std::abs(v[2]));
      return Vector3i(v / g);
    }
  }
  return std::nullopt;
}

// Index into kRotAxes that is the rotation axis (eigenvector) of a proper
// rotation; std::nullopt for the identity / when none matches. A non-identity
// proper rotation fixes exactly one line, which the table holds in one sign.
[[nodiscard]] std::optional<int> axis_index(Matrix3i const &proper_rot) {
  return primitive_axis(proper_rot).and_then([](Vector3i const &axis) {
    return table_index(axis).or_else(
        [&] { return table_index(Vector3i(-axis)); });
  });
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
  return {std::from_range,
          std::views::iota(0, kNumRotAxes) | std::views::filter([&](int i) {
            return (sum_rot * rot_axis(i)).isZero();
          })};
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

// Laue class 2/m: the two-fold axis is b; a and c are the two shortest
// distinct axes orthogonal to it.
[[nodiscard]] std::optional<AxisTriple> laue2m(PointSymmetry const &ps) {
  auto const two_fold = principal_axis(ps, 2);
  if (!two_fold) {
    return std::nullopt;
  }
  auto const ortho = orthogonal_axes(two_fold->prop_rot, 2);
  auto const first = shortest_axis(ortho);
  if (!first) {
    return std::nullopt;
  }
  auto const second = shortest_axis(
      ortho | std::views::filter([&](int idx) { return idx != *first; }));
  if (!second) {
    return std::nullopt;
  }
  return AxisTriple{{SignedAxis{*first, 1}, SignedAxis{two_fold->axis, 1},
                     SignedAxis{*second, 1}}};
}

// The in-plane axes (a, b) for a one-axis Laue class: the first orthogonal
// axis whose image under prop_rot is another, up to sign, with |det| < 4 so the
// cell is not F-centred. nullopt if none qualifies.
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

[[nodiscard]] AxisTriple sort_axes(AxisTriple axes) {
  std::ranges::sort(axes, {}, &SignedAxis::index);
  if (transformation_from_axes(axes).determinant() < 0) {
    std::swap(axes[1], axes[2]);
  }
  return axes;
}

// The aperiodic-axis component of candidate axis `idx` (rot_axes[idx][ap]).
[[nodiscard]] constexpr int aperiodic_component(int idx, int aperiodic_axis) {
  return kRotAxes[idx][aperiodic_axis];
}

// Exactly two of the three axes must lie in the periodic plane and one along
// the aperiodic axis. Move that one to c and orient for a positive
// determinant; nullopt for an inclined configuration.
[[nodiscard]] std::optional<AxisTriple> layer_sort_axes(AxisTriple axes,
                                                        int aperiodic_axis) {
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
    return std::nullopt;
  }
  int const axis_pos = *std::ranges::find_if(positions, along_aperiodic);
  std::swap(axes[static_cast<std::size_t>(axis_pos)], axes[2]);
  if (transformation_from_axes(axes).determinant() < 0) {
    std::swap(axes[0], axes[1]);
  }
  return axes;
}

// Layer LAUE2M: unlike the 3D case the two-fold axis becomes axis a (its
// position relative to the aperiodic axis distinguishes oblique vs. rectangular
// monoclinic layers); the remaining two axes are chosen accordingly.
[[nodiscard]] std::optional<AxisTriple> layer_laue2m(PointSymmetry const &ps,
                                                     int aperiodic_axis) {
  auto const two_fold = principal_axis(ps, 2);
  if (!two_fold) {
    return std::nullopt;
  }
  auto const ortho = orthogonal_axes(two_fold->prop_rot, 2);
  if (ortho.empty()) {
    return std::nullopt;
  }

  auto const pick = [&](std::optional<int> b,
                        std::optional<int> c) -> std::optional<AxisTriple> {
    if (!b || !c) {
      return std::nullopt;
    }
    return AxisTriple{
        {SignedAxis{two_fold->axis, 1}, SignedAxis{*b, 1}, SignedAxis{*c, 1}}};
  };

  int const a0 = aperiodic_component(two_fold->axis, aperiodic_axis);
  if (a0 == 1 || a0 == -1) {
    // Monoclinic/oblique: the two-fold is along the aperiodic axis; b and c
    // are the two shortest distinct orthogonal axes.
    auto const first = shortest_axis(ortho);
    if (!first) {
      return std::nullopt;
    }
    return pick(first, shortest_axis(ortho | std::views::filter([&](int idx) {
                                       return idx != *first;
                                     })));
  }
  if (a0 == 0) {
    // Monoclinic/rectangular: the second axis lies in the periodic plane, the
    // third along the aperiodic axis.
    auto const in_plane =
        shortest_axis(ortho | std::views::filter([&](int idx) {
                        return aperiodic_component(idx, aperiodic_axis) == 0;
                      }));
    auto const out_plane =
        shortest_axis(ortho | std::views::filter([&](int idx) {
                        int const c = aperiodic_component(idx, aperiodic_axis);
                        return c == 1 || c == -1;
                      }));
    return pick(in_plane, out_plane);
  }
  return std::nullopt;
}

// Laue classes mmm, m-3, m-3m. For a layer cell the three axes are reordered so
// the aperiodic axis is c.
[[nodiscard]] std::optional<AxisTriple>
lauennn(PointSymmetry const &ps, int rot_order,
        std::optional<int> aperiodic_axis) {
  // The distinct axes of the rotations of `rot_order`, in encounter order;
  // mmm / m-3 / m-3m each have exactly three.
  boost::container::static_vector<int, 3> idx;
  for (Matrix3i const &r : ps) {
    Matrix3i const prop_rot = proper_rotation(r);
    if ((prop_rot.trace() == -1 && rot_order == 2) ||
        (prop_rot.trace() == 1 && rot_order == 4)) {
      auto const axis = axis_index(prop_rot);
      if (axis && idx.size() < idx.capacity() &&
          !std::ranges::contains(idx, *axis)) {
        idx.push_back(*axis);
      }
    }
  }
  if (idx.size() != 3) {
    return std::nullopt;
  }
  AxisTriple axes;
  for (auto [a, id] : std::views::zip(axes, idx)) {
    a = {id, 1};
  }
  if (aperiodic_axis) {
    return layer_sort_axes(axes, *aperiodic_axis);
  }
  return sort_axes(axes);
}

[[nodiscard]] std::optional<AxisTriple>
get_axes(Laue laue, PointSymmetry const &ps,
         std::optional<int> aperiodic_axis) {
  switch (laue) {
  case Laue::laue_1:
    return AxisTriple{{SignedAxis{0, 1}, SignedAxis{1, 1}, SignedAxis{2, 1}}};
  case Laue::laue_2m:
    return aperiodic_axis ? layer_laue2m(ps, *aperiodic_axis) : laue2m(ps);
  case Laue::laue_mmm:
  case Laue::laue_m3:
    return lauennn(ps, 2, aperiodic_axis);
  case Laue::laue_4m:
  case Laue::laue_4mmm:
    // The 4-fold axis is the aperiodic axis for a layer; laue_one_axis already
    // places it at c, so no aperiodic-specific handling is needed.
    return laue_one_axis(ps, 4);
  case Laue::laue_3:
  case Laue::laue_3m:
  case Laue::laue_6m:
  case Laue::laue_6mmm:
    return laue_one_axis(ps, 3);
  case Laue::laue_m3m:
    return lauennn(ps, 4, aperiodic_axis);
  default:
    return std::nullopt;
  }
}

} // namespace

template <GroupFamily F>
Result<PointgroupTransform>
identify_point_group(std::span<Matrix3i const> rotations,
                     std::optional<int> layer_axis) {
  PointSymmetry const ps = unique_rotations(rotations);
  int const pg_num = pointgroup_number(ps);
  if (pg_num == 0) {
    return leaf::new_error(e_pointgroup_not_found{});
  }
  if constexpr (F == GroupFamily::layer) {
    // Layer groups have no cubic point groups (numbers 28..32).
    if (pg_num >= 28) {
      return leaf::new_error(e_pointgroup_not_found{});
    }
  } else {
    layer_axis = std::nullopt; // the 3D path never sorts for an aperiodic axis
  }
  std::optional<int> const aperiodic_axis = layer_axis;

  PointgroupTransform result;
  result.pointgroup = pointgroup_by_number(pg_num);
  auto const axes = get_axes(result.pointgroup.laue, ps, aperiodic_axis);
  if (!axes) {
    return leaf::new_error(e_pointgroup_not_found{});
  }
  result.transformation = transformation_from_axes(*axes);
  return result;
}

template Result<PointgroupTransform>
    identify_point_group<GroupFamily::space>(std::span<Matrix3i const>,
                                             std::optional<int>);
template Result<PointgroupTransform>
    identify_point_group<GroupFamily::layer>(std::span<Matrix3i const>,
                                             std::optional<int>);

} // namespace seitz::symmetry
