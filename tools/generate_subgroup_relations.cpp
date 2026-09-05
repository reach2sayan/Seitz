// Offline generator for tools/subgroup_relations.csv: re-derives the
// translationengleiche maximal-subgroup graph of the 230 space groups in-house
// (group::derive_t_subgroup_edges) and prints it as CSV. Slow -- a
// determination per candidate subgroup -- and it links this library, whose
// group/subgroup_graph.hpp includes the table transcribed from that CSV, so it
// is run by hand, never as part of the build.
//
//   cmake --build <dir> --target generate_subgroup_relations
//   ./generate_subgroup_relations > tools/subgroup_relations.csv
//
// The build then transcribes that CSV into the constexpr header
// (tools/transcribe_subgroup_relations.py).
#include <seitz/group/subgroup_graph.hpp>

#include <cstdio>

int main() {
  std::printf("super,sub,index\n");
  for (auto const &e : seitz::group::derive_t_subgroup_edges()) {
    std::printf("%d,%d,%d\n", e.super, e.sub, e.index);
  }
  return 0;
}
