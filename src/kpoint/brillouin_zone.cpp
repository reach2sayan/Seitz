#include <spglib/kpoint/brillouin_zone.hpp>

#include <spglib/kpoint/grid.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

namespace spglib::kpoint {

namespace {

// The 125 reciprocal-lattice offsets searched per grid point (the ±2 cube,
// ordered 0,1,2,-2,-1 on each axis). Kept Eigen-free (literal type) and mapped
// to Vector3i at use.
inline constexpr std::array<std::array<int, 3>, 125> kBzSearchSpace = {{
    {0, 0, 0},   {0, 0, 1},   {0, 0, 2},   {0, 0, -2},   {0, 0, -1},
    {0, 1, 0},   {0, 1, 1},   {0, 1, 2},   {0, 1, -2},   {0, 1, -1},
    {0, 2, 0},   {0, 2, 1},   {0, 2, 2},   {0, 2, -2},   {0, 2, -1},
    {0, -2, 0},  {0, -2, 1},  {0, -2, 2},  {0, -2, -2},  {0, -2, -1},
    {0, -1, 0},  {0, -1, 1},  {0, -1, 2},  {0, -1, -2},  {0, -1, -1},
    {1, 0, 0},   {1, 0, 1},   {1, 0, 2},   {1, 0, -2},   {1, 0, -1},
    {1, 1, 0},   {1, 1, 1},   {1, 1, 2},   {1, 1, -2},   {1, 1, -1},
    {1, 2, 0},   {1, 2, 1},   {1, 2, 2},   {1, 2, -2},   {1, 2, -1},
    {1, -2, 0},  {1, -2, 1},  {1, -2, 2},  {1, -2, -2},  {1, -2, -1},
    {1, -1, 0},  {1, -1, 1},  {1, -1, 2},  {1, -1, -2},  {1, -1, -1},
    {2, 0, 0},   {2, 0, 1},   {2, 0, 2},   {2, 0, -2},   {2, 0, -1},
    {2, 1, 0},   {2, 1, 1},   {2, 1, 2},   {2, 1, -2},   {2, 1, -1},
    {2, 2, 0},   {2, 2, 1},   {2, 2, 2},   {2, 2, -2},   {2, 2, -1},
    {2, -2, 0},  {2, -2, 1},  {2, -2, 2},  {2, -2, -2},  {2, -2, -1},
    {2, -1, 0},  {2, -1, 1},  {2, -1, 2},  {2, -1, -2},  {2, -1, -1},
    {-2, 0, 0},  {-2, 0, 1},  {-2, 0, 2},  {-2, 0, -2},  {-2, 0, -1},
    {-2, 1, 0},  {-2, 1, 1},  {-2, 1, 2},  {-2, 1, -2},  {-2, 1, -1},
    {-2, 2, 0},  {-2, 2, 1},  {-2, 2, 2},  {-2, 2, -2},  {-2, 2, -1},
    {-2, -2, 0}, {-2, -2, 1}, {-2, -2, 2}, {-2, -2, -2}, {-2, -2, -1},
    {-2, -1, 0}, {-2, -1, 1}, {-2, -1, 2}, {-2, -1, -2}, {-2, -1, -1},
    {-1, 0, 0},  {-1, 0, 1},  {-1, 0, 2},  {-1, 0, -2},  {-1, 0, -1},
    {-1, 1, 0},  {-1, 1, 1},  {-1, 1, 2},  {-1, 1, -2},  {-1, 1, -1},
    {-1, 2, 0},  {-1, 2, 1},  {-1, 2, 2},  {-1, 2, -2},  {-1, 2, -1},
    {-1, -2, 0}, {-1, -2, 1}, {-1, -2, 2}, {-1, -2, -2}, {-1, -2, -1},
    {-1, -1, 0}, {-1, -1, 1}, {-1, -1, 2}, {-1, -1, -2}, {-1, -1, -1},
}};

[[nodiscard]] Vector3i search_offset(std::size_t j) {
  auto const &o = kBzSearchSpace[j];
  return Vector3i(o[0], o[1], o[2]);
}

// Adaptive tolerance for merging BZ-boundary points: 0.01 * max over axes of
// |reciprocal vector|^2 / mesh^2.
[[nodiscard]] double bz_tolerance(Matrix3d const &rec_lattice,
                                  Vector3i const &mesh) {
  auto const lengths = std::views::iota(0, 3) |
                       std::views::transform([&](int i) {
                         return rec_lattice.col(i).squaredNorm() /
                                (static_cast<double>(mesh[i]) *
                                 static_cast<double>(mesh[i]));
                       });
  return 0.01 * *std::ranges::max_element(lengths);
}

} // namespace

BzGrid relocate_BZ_grid_address(std::vector<Vector3i> const &grid_address,
                                Vector3i const &mesh,
                                Matrix3d const &rec_lattice,
                                Vector3i const &is_shift) {
  Vector3i const bzmesh = mesh * 2;
  const auto num_bzmesh = static_cast<std::size_t>(bzmesh.prod());
  std::size_t const total = grid_address.size();
  double const tolerance = bz_tolerance(rec_lattice, mesh);

  BzGrid out;
  out.bz_map.assign(num_bzmesh, std::nullopt);
  // The first `total` slots hold each input point's closest image (at its
  // original index); boundary duplicates are appended afterwards.
  out.bz_grid_address.resize(total);

  for (std::size_t i = 0; i < total; ++i) {
    // Squared distance to the origin of each candidate image.
    std::array<double, kBzSearchSpace.size()> distance{};
    for (std::size_t j = 0; j < kBzSearchSpace.size(); ++j) {
      Vector3i const offset = search_offset(j);
      const auto mesh_d = mesh.cast<double>().array();

      Vector3d q = (grid_address[i].cast<double>().array() +
                    offset.cast<double>().array() * mesh_d +
                    0.5 * is_shift.cast<double>().array()) /
                   mesh_d;
      distance[j] = (rec_lattice * q).squaredNorm();
    }

    auto const min_it = std::ranges::min_element(distance);
    double const min_distance = *min_it;
    auto const min_index =
        static_cast<std::size_t>(std::ranges::distance(distance.begin(), min_it));

    for (std::size_t j = 0; j < kBzSearchSpace.size(); ++j) {
      if (distance[j] >= min_distance + tolerance) {
        continue;
      }
      const Eigen::Map<const Vector3i> shift(kBzSearchSpace[j].data());
      const Vector3i bz_address = grid_address[i] + shift.cwiseProduct(mesh);
      std::size_t gp = 0;
      if (j == min_index) {
        gp = i;
        out.bz_grid_address[i] = bz_address;
      } else {
        gp = out.bz_grid_address.size();
        out.bz_grid_address.push_back(bz_address);
      }
      Vector3i const bz_address_double = bz_address.array() * 2 + is_shift.array();
      out.bz_map[grid_point_double_mesh(bz_address_double, bzmesh)] = gp;
    }
  }

  return out;
}

std::vector<std::optional<std::size_t>> BZ_grid_points_by_rotations(
    Vector3i const &address_orig, std::vector<Matrix3i> const &rot_reciprocal,
    Vector3i const &mesh, Vector3i const &is_shift,
    std::vector<std::optional<std::size_t>> const &bz_map) {
  Vector3i const bzmesh = mesh * 2;
  Vector3i const address_double_orig =
      address_orig.array() * 2 + is_shift.array();

  std::vector<std::optional<std::size_t>> out;
  out.reserve(rot_reciprocal.size());
  std::ranges::transform(
      rot_reciprocal, std::back_inserter(out), [&](Matrix3i const &rot) {
        return bz_map[grid_point_double_mesh(
            Vector3i(rot * address_double_orig), bzmesh)];
      });
  return out;
}

} // namespace spglib::kpoint
