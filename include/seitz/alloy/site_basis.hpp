#pragma once

#include <seitz/core/types.hpp>

#include <cstddef>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::alloy {

// The per-site cluster basis: point functions theta_f(occupation) whose
// products over a cluster's points are that cluster's basis functions. A site
// admitting k species has k - 1 non-constant functions over k occupations: one
// dense (k-1 x k) block per species count.
//
// Blocks are keyed by k itself, not ATAT's "site type" = k - 2; the bias
// survives only in operator[]'s subscript arithmetic. A Chebyshev, polynomial
// or indicator basis would be another factory over this same storage.
class SiteBasis {
public:
  // van de Walle's orthonormal trigonometric basis, covering sublattices that
  // admit up to `max_species` species:
  //
  //   theta_{2t-2}(s) = -cos(2 pi s t / k),   t = 1 .. floor(k/2)
  //   theta_{2t-1}(s) = -sin(2 pi s t / k),   t = 1 .. ceil(k/2) - 1
  //
  // For k = 2 this is the Ising spin, theta_0(0) = -1 and theta_0(1) = +1.
  [[nodiscard]] static SiteBasis trigonometric(int max_species);

  // theta_function(species_index) on a site admitting `species` species.
  [[nodiscard]] double operator[](int species, int function,
                                  int occupation) const noexcept {
    return block(species)(function, occupation);
  }

  // sum_q theta_function(q)^2 over the occupations the site admits: the
  // per-point factor of the v-matrix normalization. Kept here because a row of
  // a column-major block is not contiguous, so no span can be handed out.
  [[nodiscard]] double sum_of_squares(int species, int function) const noexcept {
    return block(species).row(function).squaredNorm();
  }

  [[nodiscard]] int max_species() const noexcept {
    return static_cast<int>(blocks_.size()) + 1;
  }

private:
  SiteBasis() = default;

  [[nodiscard]] MatrixXd const &block(int species) const noexcept {
    return blocks_[static_cast<std::size_t>(species - 2)];
  }

  // blocks_[i] covers k = i + 2; a one-species site is a spectator with no
  // point function, so nothing below binary.
  std::vector<MatrixXd> blocks_;
};

} // namespace seitz::alloy

#pragma GCC visibility pop
