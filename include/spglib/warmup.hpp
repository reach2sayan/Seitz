#pragma once

#include <future>

namespace spglib {

struct WarmupOptions {
  // data::operations_from_database(hall > 0):Hit on essentially every dataset()
  bool space_group_operations = true;
  // data::operations_from_database(hall < 0): the 116 layer-group settings.
  bool layer_group_operations = true;
  // group::SubgroupGraph::instance(): builds the t-subgroup Boost graph and its
  // transitive closure.
  bool subgroup_graph = false;
};

void warmup(WarmupOptions opts = {});
[[nodiscard]] std::future<void> warmup_async(WarmupOptions opts = {});

} // namespace spglib
