#include <spglib/reduce/delaunay.hpp>

#include <spglib/math/integer_matrix.hpp>

#include <array>
#include <cmath>
#include <optional>
#include <ranges>
#include <utility>

namespace spglib::reduce {

namespace {

constexpr double kZeroPrec = 1e-10;
constexpr int kMaxAttempts = 1000; // spglib default (SPGLIB_NUM_ATTEMPTS)

// Extended Delaunay basis {b1, b2, b3, b4} with b4 = -(b1 + b2 + b3).
[[nodiscard]] std::array<Vector3d, 4> extended_basis(Matrix3d const &lattice) {
  return {lattice.col(0), lattice.col(1), lattice.col(2),
          Vector3d(-(lattice.col(0) + lattice.col(1) + lattice.col(2)))};
}

// One Delaunay reduction step on the extended basis: if any pair has a positive
// dot product, fold and flip and report "not yet reduced". Returns true when no
// positive dot product remains (delaunay.c delaunay_reduce_basis, rank 3).
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
// it has not converged within kMaxAttempts steps (delaunay.c
// delaunay_reduce_basis, rank 3).
[[nodiscard]] std::optional<std::array<Vector3d, 4>>
reduce_basis(std::array<Vector3d, 4> b, double symprec) {
  for (std::size_t attempt = 0; attempt < kMaxAttempts; ++attempt)
    if (reduce_step(b, symprec)) {
      return b;
    }
  return std::nullopt;
}

// From the candidate set {b1, b2, b1+b2, b3, b4, b2+b3, b3+b1}, pick the three
// shortest linearly independent vectors (delaunay.c
// get_delaunay_shortest_vectors, rank 3) into basis[0..2].
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

} // namespace

Result<Matrix3d> delaunay_reduce(Matrix3d const &lattice, double symprec) {
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
  return *reduced;
}

} // namespace spglib::reduce
