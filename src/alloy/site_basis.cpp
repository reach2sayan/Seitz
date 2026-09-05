#include <seitz/alloy/site_basis.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <ranges>
#include <utility>

namespace seitz::alloy {

SiteBasis SiteBasis::trigonometric(int max_species) {
  // Floored at binary: a table with no block at all would make every
  // subscript below ill-formed, and a one-species site never reaches one.
  int const largest = std::max(2, max_species);

  SiteBasis basis;
  basis.blocks_.reserve(static_cast<std::size_t>(largest - 1));
  for (int const k : std::views::iota(2, largest + 1)) {
    MatrixXd block = MatrixXd::Zero(k - 1, k);
    for (int const s : std::views::iota(0, k)) {
      // Cosines take the even rows and sines the odd ones, so the two families
      // interleave into exactly k - 1 functions however k splits.
      for (int const t : std::views::iota(1, k / 2 + 1)) {
        block(2 * t - 2, s) = -std::cos(2.0 * std::numbers::pi * s * t / k);
      }
      for (int const t : std::views::iota(1, (k + 1) / 2)) {
        block(2 * t - 1, s) = -std::sin(2.0 * std::numbers::pi * s * t / k);
      }
    }
    basis.blocks_.push_back(std::move(block));
  }
  return basis;
}

} // namespace seitz::alloy
