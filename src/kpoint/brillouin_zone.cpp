#include <cppcrystal/kpoint/brillouin_zone.hpp>

#include <cppcrystal/kpoint/grid.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ranges>
#include <vector>

namespace cppcrystal::kpoint {

namespace {

// The 125 reciprocal-lattice offsets searched per grid point: the ±2 cube in
// odometer order over the digits 0,1,2,-2,-1 per axis (x outermost, z varying
// fastest). Held as plain int triples rather than Vector3i: Eigen's (x, y, z)
// constructor is not constexpr, and its CRTP coefficient accessors are not
// constant-evaluable under MSVC, so an Eigen table cannot be built portably at
// compile time. `bz_offset` materialises a Vector3i at the use site.
[[nodiscard]] consteval std::array<std::array<int, 3>, 125> bz_search_space() {
  constexpr std::array<int, 5> digits{0, 1, 2, -2, -1};
  std::array<std::array<int, 3>, 125> table{};
  std::size_t n = 0;
  for (int const x : digits) {
    for (int const y : digits) {
      for (int const z : digits) {
        table[n++] = {x, y, z};
      }
    }
  }
  return table;
}

inline constexpr std::array<std::array<int, 3>, 125> kBzSearchSpace =
    bz_search_space();

[[nodiscard]] inline Vector3i bz_offset(std::size_t j) noexcept {
  auto const &o = kBzSearchSpace[j];
  return {o[0], o[1], o[2]};
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

  auto const mesh_d = mesh.cast<double>().array();
  auto const shift_d = 0.5 * is_shift.cast<double>().array();
  for (auto const [index, address] : grid_address | std::views::enumerate) {
    auto const i = static_cast<std::size_t>(index);
    // Squared distance to the origin of each candidate image.
    std::array<double, kBzSearchSpace.size()> distance{};
    std::ranges::transform(
        std::views::iota(std::size_t{0}, kBzSearchSpace.size()),
        distance.begin(), [&](std::size_t j) {
          Vector3d const q = (address.cast<double>().array() +
                              bz_offset(j).cast<double>().array() * mesh_d +
                              shift_d) /
                             mesh_d;
          return (rec_lattice * q).squaredNorm();
        });

    auto const min_it = std::ranges::min_element(distance);
    double const min_distance = *min_it;
    auto const min_index = static_cast<std::size_t>(min_it - distance.begin());

    // Every image within tolerance of the closest is a BZ point: the closest
    // keeps the input index, the ties are appended.
    for (auto const [offset, d] : distance | std::views::enumerate) {
      if (d >= min_distance + tolerance) {
        continue;
      }
      auto const j = static_cast<std::size_t>(offset);
      Vector3i const bz_address = address + bz_offset(j).cwiseProduct(mesh);
      std::size_t gp = i;
      if (j == min_index) {
        out.bz_grid_address[i] = bz_address;
      } else {
        gp = out.bz_grid_address.size();
        out.bz_grid_address.push_back(bz_address);
      }
      Vector3i const bz_address_double =
          bz_address.array() * 2 + is_shift.array();
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

} // namespace cppcrystal::kpoint
