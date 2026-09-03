#include <cppcrystal/kpoint/mesh.hpp>

#include <cppcrystal/core/tolerance.hpp> // approx_equal

#include "core/matrix_order.hpp" // unique_by_rotation
#include "math/fractional.hpp"   // math::round_to_int

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <vector>

namespace cppcrystal::kpoint {

namespace {

[[nodiscard]] Vector3i to_eigen(Address a) noexcept { return {a[0], a[1], a[2]}; }
[[nodiscard]] Address from_eigen(Vector3i const &v) noexcept {
  return {v[0], v[1], v[2]};
}

// True -> the fast "normal" reduction; false -> the "distortion" path (3/6-fold
// rotations or non-conventional cells). NOTE the deliberate quirk that the a=b
// and b=c flags test the same column — kept for bit-identical behaviour against
// the reference.
[[nodiscard]] bool has_conventional_symmetry(Mesh const &mesh,
                                             std::span<Matrix3i const> rots) {
  if (std::ranges::any_of(rots, [](Matrix3i const &rot) {
        return rot.cwiseAbs().sum() > 3;
      })) {
    return false;
  }

  Vector3i const e_b(0, 1, 0);
  Vector3i const e_c(0, 0, 1);
  bool const eq_ab = std::ranges::any_of( // a == b axis present
      rots, [&](Matrix3i const &rot) { return rot.col(0) == e_b; });
  bool const eq_bc = eq_ab; // b == c (intentionally the same column as a == b)
  bool const eq_ca = std::ranges::any_of( // c == a
      rots, [&](Matrix3i const &rot) { return rot.col(0) == e_c; });

  auto const &d = mesh.divisions();
  auto const &s = mesh.shift();
  return (!eq_ab || (d[0] == d[1] && s[0] == s[1])) &&
         (!eq_bc || (d[1] == d[2] && s[1] == s[2])) &&
         (!eq_ca || (d[2] == d[0] && s[2] == s[0]));
}

// Smallest grid-point index in the orbit of `doubled` under the reciprocal
// group. The group is closed, so iterating every rotation reaches the whole
// orbit and the minimum is canonical.
[[nodiscard]] std::size_t normal_representative(
    Address doubled, Mesh const &mesh, std::size_t self,
    std::span<Matrix3i const> rots) {
  Vector3i const ad = to_eigen(doubled);
  std::size_t best = self;
  for (Matrix3i const &rot : rots) {
    best = std::min(best, mesh.index_of_doubled(from_eigen(rot * ad)));
  }
  return best;
}

// As above, but for the distortion path: the rotation is applied in 64-bit
// (scaled by `divisor`) and a rotated point is only valid when it divides
// evenly and its parity matches the shift.
[[nodiscard]] std::size_t distortion_representative(
    Address doubled, Mesh const &mesh, std::size_t self,
    std::array<std::int64_t, 3> const &divisor,
    std::span<Matrix3i const> rots) {
  auto const &shift = mesh.shift();
  std::array<std::int64_t, 3> long_address{};
  for (std::size_t j = 0; j < 3; ++j) {
    long_address[j] = static_cast<std::int64_t>(doubled[j]) * divisor[j];
  }

  std::size_t best = self;
  for (auto const &rot : rots) {
    Address rotated{};
    bool indivisible = false;
    for (int k = 0; k < 3; ++k) {
      std::int64_t const value =
          static_cast<std::int64_t>(rot(k, 0)) * long_address[0] +
          static_cast<std::int64_t>(rot(k, 1)) * long_address[1] +
          static_cast<std::int64_t>(rot(k, 2)) * long_address[2];
      auto const div = divisor[static_cast<std::size_t>(k)];
      if (value % div != 0) {
        indivisible = true;
        break;
      }
      rotated[static_cast<std::size_t>(k)] = static_cast<int>(value / div);
      bool const odd = rotated[static_cast<std::size_t>(k)] % 2 != 0;
      if (odd != shift[static_cast<std::size_t>(k)]) {
        indivisible = true;
        break;
      }
    }
    if (!indivisible) {
      best = std::min(best, mesh.index_of_doubled(rotated));
    }
  }
  return best;
}

// The 125 reciprocal-lattice offsets searched per grid point: the +-2 cube in
// odometer order over the digits 0,1,2,-2,-1 per axis (x outermost, z varying
// fastest). Plain int triples rather than Vector3i: Eigen's (x, y, z)
// constructor is not constexpr, and its coefficient accessors are not
// constant-evaluable under MSVC, so an Eigen table cannot be built portably at
// compile time.
inline constexpr std::array<Address, 125> kBzSearchSpace = [] {
  constexpr std::array<int, 5> digits{0, 1, 2, -2, -1};
  std::array<Address, 125> table{};
  std::size_t n = 0;
  for (int const x : digits) {
    for (int const y : digits) {
      for (int const z : digits) {
        table[n++] = {x, y, z};
      }
    }
  }
  return table;
}();

// Adaptive tolerance for merging BZ-boundary points: 0.01 * max over axes of
// |reciprocal vector|^2 / divisions^2.
[[nodiscard]] double bz_tolerance(Matrix3d const &reciprocal,
                                  Address const &divisions) {
  auto const lengths =
      std::views::iota(0, 3) | std::views::transform([&](int i) {
        auto const d = static_cast<double>(divisions[static_cast<std::size_t>(i)]);
        return reciprocal.col(i).squaredNorm() / (d * d);
      });
  return 0.01 * *std::ranges::max_element(lengths);
}

} // namespace

ReciprocalMesh::ReciprocalMesh(Mesh mesh, std::vector<Matrix3i> rotations)
    : mesh_(mesh), rotations_(std::move(rotations)) {
  bool const normal = has_conventional_symmetry(mesh_, rotations_);
  auto const &d = mesh_.divisions();
  std::array<std::int64_t, 3> const divisor{
      static_cast<std::int64_t>(d[1]) * d[2],
      static_cast<std::int64_t>(d[2]) * d[0],
      static_cast<std::int64_t>(d[0]) * d[1],
  };

  mapping_.reserve(mesh_.size());
  for (auto const [index, address] : mesh_.addresses() | std::views::enumerate) {
    auto const i = static_cast<std::size_t>(index);
    Address const doubled = mesh_.doubled_address(address);
    mapping_.push_back(normal ? normal_representative(doubled, mesh_, i,
                                                      rotations_)
                              : distortion_representative(doubled, mesh_, i,
                                                          divisor, rotations_));
  }
  num_irreducible_ = static_cast<std::size_t>(std::ranges::count_if(
      mapping_ | std::views::enumerate,
      [](auto const &e) {
        return std::get<1>(e) == static_cast<std::size_t>(std::get<0>(e));
      }));
}

ReciprocalMesh ReciprocalMesh::from_rotations(
    Mesh mesh, std::span<Matrix3i const> real_rotations,
    TimeReversal time_reversal) {
  // All transposes first, then (with time reversal) all inversion partners
  // (-transpose); de-duplication keeps first occurrence, so order matters.
  std::vector<Matrix3i> candidates;
  candidates.reserve(real_rotations.size() *
                     (time_reversal == TimeReversal::on ? 2 : 1));
  std::ranges::transform(real_rotations, std::back_inserter(candidates),
                         [](Matrix3i const &r) { return r.transpose(); });
  if (time_reversal == TimeReversal::on) {
    std::ranges::transform(real_rotations, std::back_inserter(candidates),
                           [](Matrix3i const &r) { return -r.transpose(); });
  }
  return ReciprocalMesh{mesh, unique_by_rotation(candidates)};
}

ReciprocalMesh
ReciprocalMesh::stabilized(std::span<Vector3d const> qpoints) const {
  auto const &d = mesh_.divisions();
  double const tolerance = 0.01 / (d[0] + d[1] + d[2]);
  auto preserves_qpoints = [&](Matrix3i const &rot) {
    return std::ranges::all_of(qpoints, [&](Vector3d const &q) {
      Vector3d const q_rot = rot.cast<double>() * q;
      return std::ranges::any_of(qpoints, [&](Vector3d const &other) {
        Vector3d const diff = q_rot - other;
        return approx_equal(diff, math::round_to_int(diff).cast<double>(),
                            tolerance);
      });
    });
  };
  return ReciprocalMesh{
      mesh_, std::vector<Matrix3i>{
                 std::from_range,
                 rotations_ | std::views::filter(preserves_qpoints)}};
}

std::vector<std::size_t> ReciprocalMesh::images_of(Address address) const {
  Vector3i const doubled = to_eigen(address) * 2 + to_eigen({
      mesh_.shift()[0] ? 1 : 0, mesh_.shift()[1] ? 1 : 0,
      mesh_.shift()[2] ? 1 : 0});
  std::vector<std::size_t> out;
  out.reserve(rotations_.size());
  std::ranges::transform(rotations_, std::back_inserter(out),
                         [&](Matrix3i const &rot) {
                           return mesh_.index_of_doubled(
                               from_eigen(rot * doubled));
                         });
  return out;
}

BrillouinZone ReciprocalMesh::brillouin_zone(Lattice const &reciprocal) const {
  Matrix3d const &rec = reciprocal.matrix();
  Mesh const bzmesh = mesh_.doubled();
  auto const &d = mesh_.divisions();
  auto const &shift = mesh_.shift();
  double const tolerance = bz_tolerance(rec, d);

  std::vector<std::optional<std::size_t>> map(bzmesh.size(), std::nullopt);
  // The first mesh.size() slots hold each input point's closest image (at its
  // original index); boundary duplicates are appended afterwards.
  std::vector<Address> addresses(mesh_.size());

  for (auto const [index, address] : mesh_.addresses() | std::views::enumerate) {
    auto const i = static_cast<std::size_t>(index);
    // Squared distance to the origin of each candidate image.
    std::array<double, kBzSearchSpace.size()> distance{};
    std::ranges::transform(
        kBzSearchSpace, distance.begin(), [&](Address const &offset) {
          Vector3d q;
          for (std::size_t k = 0; k < 3; ++k) {
            auto const dk = static_cast<double>(d[k]);
            q[static_cast<int>(k)] =
                (static_cast<double>(address[k]) +
                 static_cast<double>(offset[k]) * dk + (shift[k] ? 0.5 : 0.0)) /
                dk;
          }
          return (rec * q).squaredNorm();
        });

    auto const min_it = std::ranges::min_element(distance);
    double const min_distance = *min_it;
    auto const min_index = static_cast<std::size_t>(min_it - distance.begin());

    // Every image within tolerance of the closest is a BZ point: the closest
    // keeps the input index, the ties are appended.
    for (auto const [slot, dist] : distance | std::views::enumerate) {
      if (dist >= min_distance + tolerance) {
        continue;
      }
      auto const j = static_cast<std::size_t>(slot);
      Address bz{};
      for (std::size_t k = 0; k < 3; ++k) {
        bz[k] = address[k] + kBzSearchSpace[j][k] * d[k];
      }
      std::size_t gp = i;
      if (j == min_index) {
        addresses[i] = bz;
      } else {
        gp = addresses.size();
        addresses.push_back(bz);
      }
      Address doubled{};
      for (std::size_t k = 0; k < 3; ++k) {
        doubled[k] = bz[k] * 2 + (shift[k] ? 1 : 0);
      }
      map[bzmesh.index_of_doubled(doubled)] = gp;
    }
  }

  return BrillouinZone{bzmesh, rotations_, std::move(addresses),
                       std::move(map)};
}

std::vector<std::optional<std::size_t>>
BrillouinZone::images_of(Address address) const {
  // mesh_ is already the doubled mesh.
  Vector3i const doubled = to_eigen(address) * 2 + to_eigen({
      mesh_.shift()[0] ? 1 : 0, mesh_.shift()[1] ? 1 : 0,
      mesh_.shift()[2] ? 1 : 0});
  std::vector<std::optional<std::size_t>> out;
  out.reserve(rotations_.size());
  std::ranges::transform(
      rotations_, std::back_inserter(out), [&](Matrix3i const &rot) {
        return map(mesh_.index_of_doubled(from_eigen(rot * doubled)));
      });
  return out;
}

} // namespace cppcrystal::kpoint
