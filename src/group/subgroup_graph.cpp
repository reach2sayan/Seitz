#include <seitz/group/subgroup_graph.hpp>

#include "math/integer_matrix.hpp"
#include <seitz/core/fractional.hpp>

#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/property_map/property_map.hpp>

#include <algorithm>
#include <array>
#include <optional>
#include <vector>

namespace seitz::group {

SubgroupEdge SubgroupEdge::of(int id) noexcept {
  data::SubgroupRelation const &r = detail::relation(id);
  return {.id = id,
          .super = r.super,
          .sub = detail::sub_of(r),
          .hall = *HallNumber::of(GroupFamily::space, r.hall),
          .kind = r.kind,
          .index = r.index,
          .basis = Eigen::Map<Eigen::Matrix<int, 3, 3, Eigen::RowMajor> const>(
                       r.basis.data())
                       .cast<double>() /
                   data::kTableDenominator,
          .origin = Eigen::Map<Vector3i const>(r.shift.data()).cast<double>() /
                    data::kTableDenominator};
}

std::optional<SymmetryOperation>
SubgroupEdge::in_subgroup_frame(SymmetryOperation const &op) const {
  Matrix3d const &p = basis;
  Matrix3d const p_inv = p.inverse();
  Matrix3d const r = op.rotation.cast<double>();
  Matrix3d const rotation = p_inv * r * p;
  if (!math::is_int_matrix(rotation, 1e-6)) {
    return std::nullopt;
  }
  return SymmetryOperation{
      .rotation = math::round_to_int(rotation),
      .translation = Vector3d(p_inv * (r * origin + op.translation - origin))};
}

namespace {

// A vertex the search has not reached carries no predecessor edge.
constexpr int kNoEdge = -1;

[[nodiscard]] constexpr bool in_range(int number) noexcept {
  return number >= 1 && number <= kNumSpaceGroups;
}

} // namespace

std::optional<std::vector<int>>
SubgroupGraph::path(int super, int sub, std::optional<SubgroupKind> kind) {
  if (!in_range(super) || !in_range(sub)) {
    return std::nullopt;
  }
  if (super == sub) {
    return std::vector<int>{};
  }

  // BFS over the kind view recording tree edges; the chain is read back.
  SubgroupView const view{kSubgroupGraph, detail::KindIs{kind}};
  std::array<boost::default_color_type, detail::kNumVertices> color{};
  std::array<int, detail::kNumVertices> predecessor{};
  predecessor.fill(kNoEdge);
  auto const color_map = boost::make_iterator_property_map(
      color.data(), boost::identity_property_map{});
  auto const predecessor_map = boost::make_iterator_property_map(
      predecessor.data(), boost::identity_property_map{});
  boost::breadth_first_search(
      view, super,
      boost::visitor(boost::make_bfs_visitor(boost::record_edge_predecessors(
                         predecessor_map, boost::on_tree_edge{})))
          .color_map(color_map));

  if (predecessor[static_cast<std::size_t>(sub)] == kNoEdge) {
    return std::nullopt;
  }
  std::vector<int> chain;
  for (int v = sub; v != super;) {
    int const edge = predecessor[static_cast<std::size_t>(v)];
    chain.push_back(edge);
    v = detail::source(edge, kSubgroupGraph);
  }
  std::ranges::reverse(chain);
  return chain;
}

} // namespace seitz::group
