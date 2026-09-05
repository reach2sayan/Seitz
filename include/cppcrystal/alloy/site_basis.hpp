#pragma once

#include <cppcrystal/core/types.hpp>

#include <cstddef>
#include <vector>

#pragma GCC visibility push(default)

namespace cppcrystal::alloy {

// The per-site cluster basis: the point functions theta_f(occupation) whose
// products over a cluster's points form that cluster's basis functions. A site
// admitting k species carries k - 1 non-constant point functions over k
// occupations, so the table is one dense (k-1 x k) block per species count.
//
// Rows are keyed by the SPECIES COUNT itself. ATAT keys them by "site type" =
// count - 2 and threads that biased number through every cluster record, where
// it silently changes meaning halfway through the CVM build; the bias buys
// nothing but a class of off-by-two bugs, so here it survives nowhere but the
// subscript arithmetic inside operator[].
//
// Keeping the basis behind a type is the extension point: a Chebyshev,
// polynomial or indicator basis would be a further named factory over exactly
// this storage, with nothing else in the module changing.
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

  // Sum of theta_function(q)^2 over every occupation q the site admits: the
  // per-point factor of the CVM v-matrix's normalization. Folded in here
  // because it is the only consumer that wants a whole row, and an Eigen row of
  // a column-major block is not contiguous -- handing out a span would mean
  // changing the storage order for one call site.
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

  // blocks_[i] covers species count i + 2; a site admitting one species is a
  // spectator and carries no point function, so there is no block below binary.
  std::vector<MatrixXd> blocks_;
};

} // namespace cppcrystal::alloy

#pragma GCC visibility pop
