#include <seitz/kpoint/mesh.hpp>

#include <seitz/core/tolerance.hpp> // approx_equal

#include "core/matrix_order.hpp" // unique_by_rotation
#include "core/parallel_for.hpp" // the mapping loop is the one hot parallel loop
#include "math/fractional.hpp"   // math::round_to_int

#include <boost/container/small_vector.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <vector>

namespace seitz::kpoint {

namespace {

// Grid points below which the reduction runs on the calling thread. Each point
// costs a min-reduction over up to 48 integer 3x3 products, so a few thousand
// of them is already far more work than spawning threads.
constexpr Index kMeshGrain = 4096;

[[nodiscard]] Vector3i to_eigen(Address a) noexcept {
  return {a[0], a[1], a[2]};
}
[[nodiscard]] Address from_eigen(Vector3i const &v) noexcept {
  return {v[0], v[1], v[2]};
}

// True -> fast "normal" reduction;
// false -> "distortion" path (3/6-fold rotations or non-conventional cells)
[[nodiscard]] bool has_conventional_symmetry(Mesh const &mesh,
                                             std::span<Matrix3i const> rots) {
  if (std::ranges::any_of(
          rots, [](Matrix3i const &rot) { return rot.cwiseAbs().sum() > 3; })) {
    return false;
  }

  // deliberate quirk that the a=b and b=c flags test the same column — kept for
  // bit-identical behaviour against the reference.
  Vector3i const e_b(0, 1, 0);
  Vector3i const e_c(0, 0, 1);
  bool const eq_ab = std::ranges::any_of(
      rots, [&](Matrix3i const &rot) { return rot.col(0) == e_b; });
  bool const eq_bc = eq_ab; // b == c (intentionally the same column as a == b)
  bool const eq_ca = std::ranges::any_of(
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
[[nodiscard]] std::size_t
normal_representative(Address doubled, Mesh const &mesh, std::size_t self,
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
[[nodiscard]] std::size_t
distortion_representative(Address doubled, Mesh const &mesh, std::size_t self,
                          std::array<std::int64_t, 3> const &divisor,
                          std::span<Matrix3i const> rots) {
  auto const &shift = mesh.shift();
  // Both arrays viewed in place rather than copied entry by entry.
  Eigen::Vector3<std::int64_t> const long_address =
      Eigen::Map<Eigen::Vector3i const>(doubled.data())
          .cast<std::int64_t>()
          .cwiseProduct(
              Eigen::Map<Eigen::Vector3<std::int64_t> const>(divisor.data()));

  std::size_t best = self;
  for (auto const &rot : rots) {
    // One 64-bit matrix-vector product rather than three transcribed dot
    // products. All three components are computed before the checks below,
    // where the old loop could stop at the first failing one -- same result,
    // and at most 48 rotations over 3 components it is not worth the branch.
    Eigen::Vector3<std::int64_t> const value =
        rot.cast<std::int64_t>() * long_address;

    Address rotated{};
    bool valid = true;
    for (Eigen::Index k = 0; k < 3; ++k) {
      auto const k_u = static_cast<std::size_t>(k);
      std::int64_t const div = divisor[k_u];
      if (value[k] % div != 0) {
        valid = false;
        break;
      }
      rotated[k_u] = static_cast<int>(value[k] / div);
      if ((rotated[k_u] % 2 != 0) != shift[k_u]) {
        valid = false;
        break;
      }
    }
    if (valid) {
      best = std::min(best, mesh.index_of_doubled(rotated));
    }
  }
  return best;
}

// The 125 offsets searched per grid point: the +-2 cube over 0,1,2,-2,-1 per
// axis. cartesian_product varies its last range fastest, pinned by the
// static_asserts below. Plain int triples -- Eigen fixed-size matrices are
// literal for construction and access, not arithmetic.
inline constexpr std::array<Address, 125> kBzSearchSpace = [] {
  constexpr std::array<int, 5> digits{0, 1, 2, -2, -1};
  std::array<Address, 125> table{};
  std::ranges::transform(std::views::cartesian_product(digits, digits, digits),
                         table.begin(), [](auto const &xyz) {
                           auto const [x, y, z] = xyz;
                           return Address{x, y, z};
                         });
  return table;
}();
static_assert(kBzSearchSpace.front() == Address{0, 0, 0});
static_assert(kBzSearchSpace[1] == Address{0, 0, 1});  // z varies fastest
static_assert(kBzSearchSpace[5] == Address{0, 1, 0});  // then y
static_assert(kBzSearchSpace[25] == Address{1, 0, 0}); // x is outermost
static_assert(kBzSearchSpace.back() == Address{-1, -1, -1});

// Adaptive tolerance for merging BZ-boundary points: 0.01 * max over axes of
// |reciprocal vector|^2 / divisions^2.
[[nodiscard]] double bz_tolerance(Matrix3d const &reciprocal,
                                  Address const &divisions) {
  auto const lengths =
      std::views::iota(0, 3) | std::views::transform([&](int i) {
        auto const d =
            static_cast<double>(divisions[static_cast<std::size_t>(i)]);
        return reciprocal.col(i).squaredNorm() / (d * d);
      });
  return 0.01 * *std::ranges::max_element(lengths);
}

} // namespace

ReciprocalMesh::ReciprocalMesh(Mesh mesh, std::vector<Matrix3i> rotations)
    : mesh_{mesh}, rotations_{std::move(rotations)} {
  bool const normal = has_conventional_symmetry(mesh_, rotations_);
  auto const &d = mesh_.divisions();
  std::array const divisor{
      static_cast<std::int64_t>(d[1]) * d[2],
      static_cast<std::int64_t>(d[2]) * d[0],
      static_cast<std::int64_t>(d[0]) * d[1],
  };

  // Indexed writes, not push_back: each mapping_[i] depends only on i and both
  // representative functions are pure. The largest raw iteration count in the
  // library (up to 10^6 points), and the only loop worth threading.
  mapping_.resize(mesh_.size());
  parallel_for(static_cast<Index>(mesh_.size()), kMeshGrain, [&](Index index) {
    auto const i = static_cast<std::size_t>(index);
    Address const doubled = mesh_.doubled_address(mesh_.address_of(i));
    mapping_[i] = normal ? normal_representative(doubled, mesh_, i, rotations_)
                         : distortion_representative(doubled, mesh_, i, divisor,
                                                     rotations_);
  });
  num_irreducible_ = static_cast<std::size_t>(std::ranges::count_if(
      mapping_ | std::views::enumerate, [](auto const &e) {
        const auto &[lindex, rindex] = e;
        return rindex == static_cast<std::size_t>(lindex);
      }));
}

ReciprocalMesh
ReciprocalMesh::from_rotations(Mesh mesh,
                               std::span<Matrix3i const> real_rotations,
                               TimeReversal time_reversal) {
  // All transposes first, then (with time reversal) all inversion partners
  // (-transpose); de-duplication keeps first occurrence, so order matters.
  std::vector<Matrix3i> candidates(
      std::from_range,
      real_rotations | std::views::transform([](Matrix3i const &r) {
        return Matrix3i(r.transpose());
      }));
  if (time_reversal == TimeReversal::on) {
    candidates.append_range(real_rotations |
                            std::views::transform([](Matrix3i const &r) {
                              return Matrix3i(-r.transpose());
                            }));
  }
  return ReciprocalMesh{mesh, unique_by_rotation(candidates)};
}

ReciprocalMesh
ReciprocalMesh::stabilized(std::span<Vector3d const> qpoints) const {
  auto const &d = mesh_.divisions();
  double const tolerance = 0.01 / (d[0] + d[1] + d[2]);
  auto preserves_qpoints = [&](Matrix3i const &rot) {
    Matrix3d const rot_d = rot.cast<double>(); // once per rotation, not per q
    return std::ranges::all_of(qpoints, [&](Vector3d const &q) {
      Vector3d const q_rot = rot_d * q;
      return std::ranges::any_of(qpoints, [&](Vector3d const &other) {
        Vector3d const diff = q_rot - other;
        return approx_equal(diff, math::round_to_int(diff).cast<double>(),
                            tolerance);
      });
    });
  };
  return ReciprocalMesh{
      mesh_,
      std::vector<Matrix3i>{
          std::from_range, rotations_ | std::views::filter(preserves_qpoints)}};
}

std::vector<std::size_t> ReciprocalMesh::images_of(Address address) const {
  Vector3i const doubled =
      to_eigen(address) * 2 +
      to_eigen({mesh_.shift()[0] ? 1 : 0, mesh_.shift()[1] ? 1 : 0,
                mesh_.shift()[2] ? 1 : 0});
  return {std::from_range,
          rotations_ | std::views::transform([&](Matrix3i const &rot) {
            return mesh_.index_of_doubled(from_eigen(rot * doubled));
          })};
}

BrillouinZone ReciprocalMesh::brillouin_zone(Lattice const &reciprocal) const {
  Matrix3d const &rec = reciprocal.matrix();
  Mesh const bzmesh = mesh_.doubled();
  auto const &d = mesh_.divisions();
  auto const &shift = mesh_.shift();
  double const tolerance = bz_tolerance(rec, d);

  std::vector<std::size_t> map(bzmesh.size(), BrillouinZone::kNoPoint);
  // The first mesh.size() slots hold each input point's closest image (at its
  // original index); boundary duplicates are appended afterwards.
  std::vector<Address> addresses(mesh_.size());

  // Split so the expensive half can be threaded: the distance scan is 125
  // products per grid point and depends on that point alone, while the pass
  // consuming it appends to a shared vector in order. Phase one records the
  // qualifying images and the closest; phase two replays the same sequence.
  struct Candidates {
    std::uint8_t closest = 0;
    boost::container::small_vector<std::uint8_t, 4> qualifying;
  };
  std::vector<Candidates> candidates(mesh_.size());

  parallel_for(static_cast<Index>(mesh_.size()), kMeshGrain, [&](Index index) {
    auto const i = static_cast<std::size_t>(index);
    Address const address = mesh_.address_of(i);
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
    Candidates &c = candidates[i];
    c.closest = static_cast<std::uint8_t>(min_it - distance.begin());
    for (auto const [slot, dist] : distance | std::views::enumerate) {
      if (dist < min_distance + tolerance) {
        c.qualifying.push_back(static_cast<std::uint8_t>(slot));
      }
    }
  });

  // Every image within tolerance of the closest is a BZ point: the closest
  // keeps the input index, the ties are appended.
  for (std::size_t i = 0; i < mesh_.size(); ++i) {
    Address const address = mesh_.address_of(i);
    Candidates const &c = candidates[i];
    for (std::uint8_t const slot : c.qualifying) {
      auto const j = static_cast<std::size_t>(slot);
      Address bz{};
      for (std::size_t k = 0; k < 3; ++k) {
        bz[k] = address[k] + kBzSearchSpace[j][k] * d[k];
      }
      std::size_t gp = i;
      if (slot == c.closest) {
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
  Vector3i const doubled =
      to_eigen(address) * 2 +
      to_eigen({mesh_.shift()[0] ? 1 : 0, mesh_.shift()[1] ? 1 : 0,
                mesh_.shift()[2] ? 1 : 0});
  return {std::from_range,
          rotations_ | std::views::transform([&](Matrix3i const &rot) {
            return map(mesh_.index_of_doubled(from_eigen(rot * doubled)));
          })};
}

} // namespace seitz::kpoint
