#include <spglib/reduce/delaunay.hpp>

#include <spglib/math/integer_matrix.hpp>

#include <array>
#include <cmath>
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
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 4; ++j)
      if (b[static_cast<std::size_t>(i)].dot(b[static_cast<std::size_t>(j)]) >
          symprec) {
        for (int k = 0; k < 4; ++k) {
          if (k != i && k != j) {
            b[static_cast<std::size_t>(k)] += b[static_cast<std::size_t>(i)];
          }
        }
        b[static_cast<std::size_t>(i)] = -b[static_cast<std::size_t>(i)];
        return false;
      }
  return true;
}

// From the candidate set {b1, b2, b1+b2, b3, b4, b2+b3, b3+b1}, pick the three
// shortest linearly independent vectors (delaunay.c
// get_delaunay_shortest_vectors, rank 3) and store them back into basis[0..2].
void shortest_vectors(std::array<Vector3d, 4> &basis, double symprec) {
  std::array<Vector3d, 7> b{basis[0],           basis[1], basis[0] + basis[1],
                            basis[2],           basis[3], basis[1] + basis[2],
                            basis[2] + basis[0]};
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      if (b[static_cast<std::size_t>(j)].squaredNorm() >
          b[static_cast<std::size_t>(j + 1)].squaredNorm() + kZeroPrec) {
        std::swap(b[static_cast<std::size_t>(j)],
                  b[static_cast<std::size_t>(j + 1)]);
      }
    }
  }

  for (int i = 2; i < 7; ++i) {
    Matrix3d m;
    m.col(0) = b[0];
    m.col(1) = b[1];
    m.col(2) = b[static_cast<std::size_t>(i)];
    if (std::abs(m.determinant()) > symprec) {
      basis[2] = b[static_cast<std::size_t>(i)];
      basis[0] = b[0];
      basis[1] = b[1];
      return;
    }
  }
}

} // namespace

Result<Matrix3d> delaunay_reduce(Matrix3d const &lattice, double symprec) {
  auto basis = extended_basis(lattice);

  bool reduced = false;
  for (int attempt = 0; attempt < kMaxAttempts; ++attempt)
    if (reduce_step(basis, symprec)) {
      reduced = true;
      break;
    }
  if (!reduced) {
    return leaf::new_error(e_delaunay_failed{});
  }
  shortest_vectors(basis, symprec);

  Matrix3d red;
  red.col(0) = basis[0];
  red.col(1) = basis[1];
  red.col(2) = basis[2];

  double const volume = red.determinant();
  if (std::abs(volume) < symprec)
    return leaf::new_error(e_delaunay_failed{});
  if (volume < 0.0)
    red = -red;

  // The change of basis from the input to the reduced lattice must be
  // unimodular.
  auto const red_inv = math::inverse(red, symprec);
  if (!red_inv) {
    return leaf::new_error(e_delaunay_failed{});
  }
  Matrix3i const change = math::round_to_int(Matrix3d(*red_inv * lattice));
  if (std::abs(change.determinant()) != 1) {
    return leaf::new_error(e_delaunay_failed{});
  }
  return red;
}

} // namespace spglib::reduce
