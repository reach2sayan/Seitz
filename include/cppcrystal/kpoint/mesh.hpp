#pragma once

#include <cppcrystal/core/lattice.hpp>
#include <cppcrystal/core/symmetry_operation.hpp>
#include <cppcrystal/core/types.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

// Reciprocal-space sampling, 3D path. A Mesh is pure grid geometry — the
// mapping between integer addresses and linear grid-point indices — and is
// constexpr on std::array<int, 3>, with Eigen only at the boundary. The
// double-mesh convention is q = (address * 2 + shift) / (divisions * 2).
namespace cppcrystal::kpoint {

using Address = std::array<int, 3>;

class BrillouinZone;

class Mesh {
public:
  // A mesh must be strictly positive on every axis, or the modulo and index
  // arithmetic below is undefined. Rejecting it here is what used to be
  // e_invalid_mesh raised at each root.
  [[nodiscard]] static constexpr std::optional<Mesh>
  of(Address divisions, std::array<bool, 3> shift = {}) noexcept {
    if (divisions[0] <= 0 || divisions[1] <= 0 || divisions[2] <= 0) {
      return std::nullopt;
    }
    return Mesh{divisions, shift};
  }

  [[nodiscard]] constexpr Address const &divisions() const noexcept {
    return divisions_;
  }
  [[nodiscard]] constexpr std::array<bool, 3> const &shift() const noexcept {
    return shift_;
  }

  [[nodiscard]] constexpr std::size_t size() const noexcept {
    return static_cast<std::size_t>(divisions_[0]) *
           static_cast<std::size_t>(divisions_[1]) *
           static_cast<std::size_t>(divisions_[2]);
  }

  // Linear index of an address, folded into [0, divisions) first:
  //   index = a0 + a1*d0 + a2*d0*d1.
  [[nodiscard]] constexpr std::size_t index_of(Address address) const noexcept {
    std::size_t stride = 1;
    std::size_t index = 0;
    for (std::size_t i = 0; i < 3; ++i) {
      index += static_cast<std::size_t>(modulo(address[i], divisions_[i])) *
               stride;
      stride *= static_cast<std::size_t>(divisions_[i]);
    }
    return index;
  }

  // Linear index of a doubled-mesh address: halve (rounding down), then fold.
  [[nodiscard]] constexpr std::size_t
  index_of_doubled(Address doubled) const noexcept {
    return index_of({floor_half(doubled[0]), floor_half(doubled[1]),
                     floor_half(doubled[2])});
  }

  // The address at a linear index, folded onto the parallelepiped.
  [[nodiscard]] constexpr Address address_of(std::size_t index) const noexcept {
    Address out{};
    for (std::size_t i = 0; i < 3; ++i) {
      auto const d = static_cast<std::size_t>(divisions_[i]);
      out[i] = static_cast<int>(index % d);
      index /= d;
      // Fold onto the parallelepiped: past the half-way point is the negative
      // image.
      if (out[i] > divisions_[i] / 2) {
        out[i] -= divisions_[i];
      }
    }
    return out;
  }

  // Every address, in grid-point-index order.
  [[nodiscard]] auto addresses() const {
    return std::views::iota(std::size_t{0}, size()) |
           std::views::transform(
               [this](std::size_t i) { return address_of(i); });
  }

  // The doubled-mesh address of a single-mesh one: 2*address + shift, folded
  // onto the doubled parallelepiped.
  [[nodiscard]] constexpr Address
  doubled_address(Address address) const noexcept {
    Address out{};
    for (std::size_t i = 0; i < 3; ++i) {
      out[i] = 2 * address[i] + (shift_[i] ? 1 : 0);
      if (out[i] > divisions_[i]) {
        out[i] -= 2 * divisions_[i];
      }
    }
    return out;
  }

  // The mesh at twice the resolution, which the Brillouin-zone map is indexed
  // on.
  [[nodiscard]] constexpr Mesh doubled() const noexcept {
    return Mesh{{2 * divisions_[0], 2 * divisions_[1], 2 * divisions_[2]},
                shift_};
  }

  [[nodiscard]] friend constexpr bool operator==(Mesh const &,
                                                 Mesh const &) = default;

private:
  constexpr Mesh(Address divisions, std::array<bool, 3> shift) noexcept
      : divisions_(divisions), shift_(shift) {}

  // Euclidean modulo into [0, m), and floor division by 2 (which differs from
  // C++ truncation for negative addresses).
  [[nodiscard]] static constexpr int modulo(int v, int m) noexcept {
    int const r = v % m;
    return r < 0 ? r + m : r;
  }
  [[nodiscard]] static constexpr int floor_half(int v) noexcept {
    return v >= 0 ? v / 2 : -((-v + 1) / 2);
  }

  Address divisions_{1, 1, 1};
  std::array<bool, 3> shift_{};
};

// Unit checks that used to need a test run.
static_assert(!Mesh::of({0, 4, 4}).has_value());
static_assert(!Mesh::of({4, -2, 4}).has_value());
static_assert(Mesh::of({4, 4, 4})->size() == 64);
static_assert(Mesh::of({4, 4, 4})->index_of({1, 2, 3}) == 1 + 2 * 4 + 3 * 16);
static_assert(Mesh::of({4, 4, 4})->index_of({-1, 0, 0}) == 3); // folds
static_assert(Mesh::of({4, 4, 4})->address_of(57) == Address{1, 2, -1});
static_assert(Mesh::of({4, 4, 4})->index_of(
                  Mesh::of({4, 4, 4})->address_of(57)) == 57); // round-trips
static_assert(Mesh::of({4, 4, 4})->doubled().divisions() == Address{8, 8, 8});
static_assert(Mesh::of({4, 4, 4}, {true, false, false})
                  ->doubled_address({1, 0, 0}) == Address{3, 0, 0});

// The symmetry reduction of a sampling mesh: the reciprocal point group and,
// for each grid point, the index of its irreducible representative (the
// smallest index in its orbit). Reciprocal rotations are the transpose of the
// real-space ones; with time reversal the inversion partner is added too.
class ReciprocalMesh {
public:
  [[nodiscard]] static ReciprocalMesh
  from_rotations(Mesh mesh, std::span<Matrix3i const> real_rotations,
                 TimeReversal time_reversal);

  // Reduced only by the rotations that map the q-point set onto itself.
  [[nodiscard]] ReciprocalMesh
  stabilized(std::span<Vector3d const> qpoints) const;

  [[nodiscard]] Mesh const &mesh() const noexcept { return mesh_; }
  [[nodiscard]] std::span<Matrix3i const> rotations() const noexcept {
    return rotations_;
  }

  // mapping()[i] == i  <=>  i is an irreducible representative.
  [[nodiscard]] std::span<std::size_t const> mapping() const noexcept {
    return mapping_;
  }
  [[nodiscard]] std::size_t num_irreducible() const noexcept {
    return num_irreducible_;
  }

  // The grid-point indices one address maps to under each rotation.
  [[nodiscard]] std::vector<std::size_t> images_of(Address address) const;

  // Relocate the grid into the first Brillouin zone of `reciprocal` (columns =
  // reciprocal basis vectors).
  [[nodiscard]] BrillouinZone
  brillouin_zone(Lattice const &reciprocal) const;

private:
  ReciprocalMesh(Mesh mesh, std::vector<Matrix3i> rotations);

  Mesh mesh_;
  std::vector<Matrix3i> rotations_;
  std::vector<std::size_t> mapping_;
  std::size_t num_irreducible_ = 0;
};

// A sampling grid relocated into the first Brillouin zone: each point moved to
// the reciprocal-lattice image closest to the origin, with boundary points
// (several equidistant images within tolerance) duplicated. The map is indexed
// on the doubled mesh.
class BrillouinZone {
public:
  // The in-BZ grid points: the first mesh.size() are the closest image of each
  // input point at its original index, followed by the boundary duplicates.
  [[nodiscard]] std::span<Address const> addresses() const noexcept {
    return addresses_;
  }

  // Doubled-mesh grid-point index -> index into addresses(); std::nullopt
  // where no BZ point falls.
  [[nodiscard]] std::optional<std::size_t>
  map(std::size_t doubled_index) const noexcept {
    return doubled_index < map_.size() ? map_[doubled_index] : std::nullopt;
  }

  // The BZ grid-point indices one address maps to under each rotation: rotate
  // on the doubled mesh, then look the result up in the map.
  [[nodiscard]] std::vector<std::optional<std::size_t>>
  images_of(Address address) const;

private:
  friend class ReciprocalMesh;
  BrillouinZone(Mesh mesh, std::span<Matrix3i const> rotations,
                std::vector<Address> addresses,
                std::vector<std::optional<std::size_t>> map)
      : mesh_(mesh), rotations_(rotations.begin(), rotations.end()),
        addresses_(std::move(addresses)), map_(std::move(map)) {}

  Mesh mesh_;
  std::vector<Matrix3i> rotations_;
  std::vector<Address> addresses_;
  std::vector<std::optional<std::size_t>> map_;
};

} // namespace cppcrystal::kpoint
