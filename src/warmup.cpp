#include <cppcrystal/warmup.hpp>

#include <cppcrystal/data/spg_database.hpp>
#include <cppcrystal/group/subgroup_graph.hpp>

#include <future>
#include <vector>

namespace cppcrystal {

namespace {

// Touch each requested singleton on its own thread so independent primers run
// concurrently; the touch discards the returned reference — the side effect is
// the one-time construction of the function-local static.
void prime(WarmupOptions opts) {
  std::vector<std::future<void>> primers;
  if (opts.space_group_operations) {
    primers.push_back(std::async(
        std::launch::async, [] { (void)data::operations_from_database(1); }));
  }
  if (opts.layer_group_operations) {
    primers.push_back(std::async(
        std::launch::async, [] { (void)data::operations_from_database(-1); }));
  }
  if (opts.subgroup_graph) {
    primers.push_back(std::async(
        std::launch::async, [] { (void)group::SubgroupGraph::instance(); }));
  }
  for (auto &f : primers) {
    f.get();
  }
}

} // namespace

void warmup(WarmupOptions opts) { prime(opts); }
std::future<void> warmup_async(WarmupOptions opts) {
  return std::async(std::launch::async, [opts] { prime(opts); });
}

} // namespace cppcrystal
