#pragma once

#include <cppcrystal/core/types.hpp>

#include <algorithm>
#include <cstddef>
#include <thread>
#include <vector>

namespace cppcrystal {

// Run `body(i)` for every i in [0, n), split into one contiguous chunk per
// hardware thread.
//
// std::jthread rather than <execution>: libstdc++ resolves par/par_unseq
// through a TBB backend selected by `__has_include(<tbb/tbb.h>)`, so the same
// source silently runs parallel on a machine with TBB headers and serial on one
// without -- and cmake/Dependencies.cmake states that every dependency is
// fetched rather than taken from the system. Threads::Threads is already a
// PUBLIC link of the library.
//
// Contract on `body`:
//   - safe to call concurrently for distinct i (callers here write only to
//     their own index, so the result is independent of the thread count and of
//     where the chunks fall);
//   - does not throw -- an exception escaping a worker terminates.
//
// Below `grain` total iterations the loop runs on the calling thread: spawning
// threads costs more than a small loop, and that is also the single-core path.
template <class F> void parallel_for(Index n, Index grain, F const &body) {
  if (n <= 0) {
    return;
  }
  auto const hardware =
      static_cast<Index>(std::max(1U, std::thread::hardware_concurrency()));
  if (n < grain || hardware <= 1) {
    for (Index i = 0; i < n; ++i) {
      body(i);
    }
    return;
  }

  Index const per_chunk = (n + hardware - 1) / hardware;
  std::vector<std::jthread> workers;
  workers.reserve(static_cast<std::size_t>(hardware) - 1);
  for (Index lo = per_chunk; lo < n; lo += per_chunk) {
    Index const hi = std::min(n, lo + per_chunk);
    workers.emplace_back([&body, lo, hi] {
      for (Index i = lo; i < hi; ++i) {
        body(i);
      }
    });
  }

  // The calling thread takes the first chunk rather than idling on the join.
  for (Index i = 0, hi = std::min(n, per_chunk); i < hi; ++i) {
    body(i);
  }
  // ~jthread joins.
}

} // namespace cppcrystal
