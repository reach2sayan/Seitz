#pragma once

#include <algorithm>
#include <cstddef>
#include <generator>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

// Private to src/alloy: the two combinatorial walks the cluster and CVM
// builders share. Coroutines rather than the callback-taking recursive lambdas
// the reference implementation uses -- three of those (one per call site,
// duplicated across two translation units) collapse into these two, and the
// call sites become ordinary range-for loops.
//
// Both yield a span into the coroutine's own buffer: valid until the generator
// is resumed, which is exactly how every caller here consumes them.
namespace seitz::alloy::detail {

// Every tuple of [0, radix[0]) x ... x [0, radix[k-1]), digit 0 varying
// FASTEST.
//
// The carry is written out rather than expressed as an algorithm: it IS the
// odometer, and the standard library has no counterpart to borrow. The order is
// load-bearing, not incidental -- it decides which member of a symmetry orbit
// is generated first and therefore becomes that orbit's stored representative,
// and it fixes the order of the correlation sums. An empty radix yields exactly
// one empty tuple, which is what makes the empty cluster fall out of the same
// loop as every other.
[[nodiscard]] inline std::generator<std::span<int const>>
mixed_radix(std::span<int const> radix) {
  if (std::ranges::any_of(radix, [](int bound) { return bound <= 0; })) {
    co_return;
  }
  std::vector<int> digits(radix.size(), 0);
  for (;;) {
    co_yield std::span<int const>{digits};
    std::size_t digit = 0;
    for (; digit < digits.size(); ++digit) {
      if (++digits[digit] < radix[digit]) {
        break;
      }
      digits[digit] = 0;
    }
    if (digit == digits.size()) {
      co_return;
    }
  }
}

// Every `choose`-subset of [0, n), in lexicographic order.
//
// The classic mask walk: a selection mask of `choose` leading ones, stepped
// with std::prev_permutation, visits every subset exactly once and hands them
// over in lexicographic order of the indices they select. That is the whole
// algorithm -- no hand-rolled "advance the rightmost index that is not at its
// ceiling" bookkeeping, which is the part that is easy to get subtly wrong.
// `choose == 0` yields one empty subset; an impossible request yields none.
[[nodiscard]] inline std::generator<std::span<std::size_t const>>
combinations(std::size_t n, int choose) {
  if (choose < 0 || std::cmp_greater(choose, n)) {
    co_return;
  }
  auto const width = static_cast<std::size_t>(choose);
  std::vector<char> mask(n, 0);
  std::ranges::fill_n(mask.begin(), static_cast<std::ptrdiff_t>(width), 1);

  std::vector<std::size_t> chosen;
  chosen.reserve(width);
  do {
    auto selected = std::views::iota(std::size_t{0}, n) |
                    std::views::filter([&mask](std::size_t i) {
                      return mask[i] != 0;
                    });
    chosen.assign(selected.begin(), selected.end());
    co_yield std::span<std::size_t const>{chosen};
  } while (std::prev_permutation(mask.begin(), mask.end()));
}

} // namespace seitz::alloy::detail
