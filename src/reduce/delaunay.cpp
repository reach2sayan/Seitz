#include <cppcrystal/core/lattice.hpp>

#include <cppcrystal/core/tolerance.hpp>
#include "math/integer_matrix.hpp"

#include <array>
#include <cmath>
#include <optional>
#include <ranges>
#include <utility>

namespace cppcrystal {

namespace {

constexpr int kMaxAttempts = 1000;

// Extended Delaunay basis {b1, b2, b3, b4} with b4 = -(b1 + b2 + b3).
[[nodiscard]] std::array<Vector3d, 4> extended_basis(Matrix3d const &lattice) {
  return {lattice.col(0), lattice.col(1), lattice.col(2),
          Vector3d(-(lattice.col(0) + lattice.col(1) + lattice.col(2)))};
}

// One Delaunay reduction step on the extended basis: if any pair has a positive
// dot product, fold and flip and report "not yet reduced". Returns true when no
// positive dot product remains.
[[nodiscard]] bool reduce_step(std::array<Vector3d, 4> &b, double symprec) {
  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = i + 1; j < 4; ++j)
      if (b[i].dot(b[j]) > symprec) {
        for (std::size_t k = 0; k < 4; ++k) {
          if (k != i && k != j) {
            b[k] += b[i];
          }
        }
        b[i] = -b[i];
        return false;
      }
  return true;
}

// Fold the extended basis until no pair has a positive dot product; nullopt if
// it has not converged within kMaxAttempts steps.
[[nodiscard]] std::optional<std::array<Vector3d, 4>>
reduce_basis(std::array<Vector3d, 4> b, double symprec) {
  for (std::size_t attempt = 0; attempt < kMaxAttempts; ++attempt)
    if (reduce_step(b, symprec)) {
      return b;
    }
  return std::nullopt;
}

// From the candidate set {b1, b2, b1+b2, b3, b4, b2+b3, b3+b1}, pick the three
// shortest linearly independent vectors into basis[0..2].
[[nodiscard]] std::array<Vector3d, 4>
shortest_vectors(std::array<Vector3d, 4> basis, double symprec) {
  std::array<Vector3d, 7> b{basis[0],           basis[1], basis[0] + basis[1],
                            basis[2],           basis[3], basis[1] + basis[2],
                            basis[2] + basis[0]};
  std::ranges::sort(b, {}, [](const Vector3d &v) { return v.squaredNorm(); });

  const auto it =
      std::ranges::find_if(std::views::drop(b, 2), [&](const Vector3d &v) {
        Matrix3d m;
        m.col(0) = b[0];
        m.col(1) = b[1];
        m.col(2) = v;
        return std::abs(m.determinant()) > symprec;
      });

  if (it != b.end()) {
    basis[0] = b[0];
    basis[1] = b[1];
    basis[2] = *it;
  }
  return basis;
}

// Assemble the reduced lattice from the first three reduced basis vectors,
// flipped to be right-handed; nullopt if the cell is degenerate.
[[nodiscard]] std::optional<Matrix3d>
oriented_cell(std::array<Vector3d, 4> const &basis, double symprec) {
  Matrix3d red;
  red << basis[0], basis[1], basis[2];

  double const volume = red.determinant();
  if (std::abs(volume) < symprec) {
    return std::nullopt;
  }
  return volume < 0.0 ? Matrix3d(-red) : red;
}

// Pass the reduced lattice through iff the change of basis from the input is
// unimodular; nullopt otherwise.
[[nodiscard]] std::optional<Matrix3d>
if_unimodular(Matrix3d const &red, Matrix3d const &lattice, double symprec) {
  return math::inverse(red, symprec)
      .and_then([&](Matrix3d const &red_inv) -> std::optional<Matrix3d> {
        Matrix3i const change = math::round_to_int(Matrix3d(red_inv * lattice));
        if (std::abs(change.determinant()) != 1) {
          return std::nullopt;
        }
        return red;
      });
}

// ---- 2D Delaunay (space-group monoclinic path, lattice rank 2) ----

// Extended 2D basis {p, q, -(p+q)} for the plane spanned by columns j, k.
// Overload of extended_basis distinguished by its two-vector argument list.
[[nodiscard]] std::array<Vector3d, 3> extended_basis(Vector3d const &p,
                                                     Vector3d const &q) {
  return {p, q, Vector3d(-(p + q))};
}

// One rank-2 Delaunay step: if any pair has a positive dot product, fold the
// third vector by +2*b[i], flip b[i], and report "not yet reduced". The rank-2
// / rank-3 overloads are distinguished by the basis array size.
[[nodiscard]] bool reduce_step(std::array<Vector3d, 3> &b, double symprec) {
  for (std::size_t i = 0; i < 2; ++i)
    for (std::size_t j = i + 1; j < 3; ++j)
      if (b[i].dot(b[j]) > symprec) {
        for (std::size_t k = 0; k < 3; ++k)
          if (k != i && k != j) {
            b[k] += 2 * b[i];
            break;
          }
        b[i] = -b[i];
        return false;
      }
  return true;
}

[[nodiscard]] std::optional<std::array<Vector3d, 3>>
reduce_basis(std::array<Vector3d, 3> b, double symprec) {
  for (std::size_t attempt = 0; attempt < kMaxAttempts; ++attempt)
    if (reduce_step(b, symprec)) {
      return b;
    }
  return std::nullopt;
}

// From {b0, b1, b2, b0+b1} pick the two shortest in-plane vectors that, with
// `unique_vec`, span a non-degenerate cell. Returns {first, second}. Overload
// distinguished from the rank-3 shortest_vectors by the `unique_vec` argument.
[[nodiscard]] std::array<Vector3d, 2>
shortest_vectors(std::array<Vector3d, 3> const &basis,
                 Vector3d const &unique_vec, double symprec) {
  std::array<Vector3d, 4> b{basis[0], basis[1], basis[2],
                            Vector3d(basis[0] + basis[1])};

  // Bubble sort ascending by squared norm, preserving order on near-ties: a
  // plain std::ranges::sort would reorder near-equal vectors and change which
  // one is picked below.
  for (std::size_t pass = 0; pass < 3; ++pass)
    for (std::size_t j = 0; j < 3; ++j)
      if (sqnorm_longer(b[j].squaredNorm(), b[j + 1].squaredNorm())) {
        std::swap(b[j], b[j + 1]);
      }

  std::array<Vector3d, 2> out{b[0], b[0]};
  for (std::size_t i = 1; i < 4; ++i) {
    Matrix3d m;
    m.col(0) = b[0];
    m.col(1) = unique_vec;
    m.col(2) = b[i];
    if (std::abs(m.determinant()) > symprec) {
      out[1] = b[i];
      break;
    }
  }
  return out;
}

} // namespace

Result<Lattice> Lattice::delaunay_in_plane(int unique_axis, double symprec) const {
  Matrix3d const &lattice = basis_;
  // The two in-plane axes (j < k) are those other than the unique axis.
  static constexpr std::array<std::array<int, 2>, 3> planes{{
      {{1, 2}},
      {{0, 2}},
      {{0, 1}},
  }};

  const auto [j, k] = planes[static_cast<std::size_t>(unique_axis)];
  Vector3d const unique_vec = lattice.col(unique_axis);

  auto const reduced =
      reduce_basis(extended_basis(lattice.col(j), lattice.col(k)), symprec)
          .transform([&](std::array<Vector3d, 3> const &b) {
            return shortest_vectors(b, unique_vec, symprec);
          });
  if (!reduced) {
    return leaf::new_error(e_delaunay_failed{});
  }

  Matrix3d red;
  red.col(unique_axis) = unique_vec;
  red.col(j) = (*reduced)[0];
  red.col(k) = (*reduced)[1];

  double const volume = red.determinant();
  if (std::abs(volume) < symprec) {
    return leaf::new_error(e_delaunay_failed{});
  }
  if (volume < 0.0) {
    red.col(unique_axis) = -red.col(unique_axis);
  }
  return Lattice{red};
}

Result<Lattice> Lattice::delaunay(double symprec) const {
  Matrix3d const &lattice = basis_;
  auto const reduced = reduce_basis(extended_basis(lattice), symprec)
                           .transform([&](std::array<Vector3d, 4> const &b) {
                             return shortest_vectors(b, symprec);
                           })
                           .and_then([&](std::array<Vector3d, 4> const &b) {
                             return oriented_cell(b, symprec);
                           })
                           .and_then([&](Matrix3d const &red) {
                             return if_unimodular(red, lattice, symprec);
                           });

  if (!reduced) {
    return leaf::new_error(e_delaunay_failed{});
  }
  return Lattice{*reduced};
}

} // namespace cppcrystal
