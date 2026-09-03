#pragma once

#include <future>

namespace cppcrystal {

struct WarmupOptions {
  // data::operations_from_database, 3D: hit on essentially every dataset().
  bool space_group_operations = true;
  // data::operations_from_database, the 116 layer-group settings.
  bool layer_group_operations = true;
};

void warmup(WarmupOptions opts = {});
[[nodiscard]] std::future<void> warmup_async(WarmupOptions opts = {});

} // namespace cppcrystal
