#include <cppcrystal/warmup.hpp>

#include <cppcrystal/core/keys.hpp>
#include <cppcrystal/group/space_group.hpp>

#include <future>
#include <ranges>
#include <vector>

namespace cppcrystal {

namespace {

// Touch every setting of one family. The touch discards the returned
// reference — the side effect is the one-time construction of the flyweight.
void prime_family(GroupFamily family) {
  for (int index : std::views::iota(1, hall_settings(family) + 1)) {
    (void)group::SpaceGroup::of(*HallNumber::of(family, index));
  }
}

void prime(Warm what) {
  std::vector<std::future<void>> primers;
  if (contains(what, Warm::space_groups)) {
    primers.push_back(std::async(std::launch::async,
                                 [] { prime_family(GroupFamily::space); }));
  }
  if (contains(what, Warm::layer_groups)) {
    primers.push_back(std::async(std::launch::async,
                                 [] { prime_family(GroupFamily::layer); }));
  }
  for (auto &f : primers) {
    f.get();
  }
}

} // namespace

void warmup(Warm what) { prime(what); }

std::future<void> warmup_async(Warm what) {
  return std::async(std::launch::async, [what] { prime(what); });
}

} // namespace cppcrystal
