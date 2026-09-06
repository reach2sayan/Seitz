#pragma once

#include <seitz/core/keys.hpp>
#include <seitz/core/symmetry_operation.hpp>
#include <seitz/core/types.hpp>
#include <seitz/data/detail/lookup.hpp>
#include <seitz/data/spg_database.hpp>
#include <seitz/data/subgroup_relations.hpp>

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/property_map/function_property_map.hpp>
#include <boost/property_map/property_map.hpp>

#include <Eigen/Core>

#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#pragma GCC visibility push(default)

namespace seitz::group {

inline constexpr int kNumSpaceGroups = 230;
using data::SubgroupKind;

// H = `sub` (setting `hall`) maximal in G = `super` (first Hall setting):
//   (a_H b_H c_H) = (a_G b_G c_G) basis,   x_G = basis x_H + origin,
// basis rational (|det| = 1/m_G for a centred parent's primitive subcell),
// origin defined modulo basis Z^3.
struct SubgroupEdge {
  int id;
  int super;
  int sub;
  HallNumber hall;
  SubgroupKind kind;
  int index;
  Matrix3d basis;
  Vector3d origin;

  [[nodiscard]] static SubgroupEdge of(int id) noexcept;

  // (R, t) of the G frame as (P^-1 R P, P^-1 (R s + t - s)); nullopt when
  // P^-1 R P is not integral (R does not preserve T_H).
  [[nodiscard]] std::optional<SymmetryOperation>
  in_subgroup_frame(SymmetryOperation const &op) const;
};

namespace detail {

inline constexpr std::size_t kNumEdges = data::kSubgroupRelations.size();
// Vertex v is space-group number v; vertex 0 is isolated.
inline constexpr std::size_t kNumVertices = kNumSpaceGroups + 1;

[[nodiscard]] constexpr data::SubgroupRelation const &
relation(int edge) noexcept {
  return data::kSubgroupRelations[static_cast<std::size_t>(edge)];
}

// The subgroup's number, off its Hall setting.
[[nodiscard]] constexpr int sub_of(data::SubgroupRelation const &r) noexcept {
  return data::spacegroup_type(*HallNumber::of(GroupFamily::space, r.hall))
      .number;
}

// Every stored Hall index names a setting (the `*` above is total).
static_assert(std::ranges::all_of(data::kSubgroupRelations, [](auto const &r) {
  return HallNumber::of(GroupFamily::space, r.hall).has_value();
}));

// CSR over the edge table: edge ids counting-sorted by super (out) and by sub
// (in) at compile time; the table is the one storage.
using EdgeIndex = data::detail::BucketIndex<int, kNumEdges, kNumVertices>;

template <int (*Endpoint)(data::SubgroupRelation const &)>
[[nodiscard]] consteval EdgeIndex group_by() {
  return data::detail::bucket_index<kNumEdges, kNumVertices, int>(
      [](int id) { return Endpoint(relation(id - 1)); },
      [](int id) { return id - 1; });
}

[[nodiscard]] constexpr int super_of(data::SubgroupRelation const &r) noexcept {
  return r.super;
}

inline constexpr EdgeIndex kOutEdges = group_by<&super_of>();
inline constexpr EdgeIndex kInEdges = group_by<&sub_of>();

using ReachabilityRow = std::bitset<kNumVertices>;

// Warshall closure of the relation (self-loops harmless): is_subgroup is a bit.
[[nodiscard]] consteval std::array<ReachabilityRow, kNumVertices>
reachability() {
  std::array<ReachabilityRow, kNumVertices> reach{};
  for (auto const &rel : data::kSubgroupRelations) {
    reach[static_cast<std::size_t>(rel.super)].set(
        static_cast<std::size_t>(sub_of(rel)));
  }
  for (std::size_t k = 1; k < kNumVertices; ++k) {
    for (std::size_t i = 1; i < kNumVertices; ++i) {
      if (reach[i].test(k)) {
        reach[i] |= reach[k];
      }
    }
  }
  return reach;
}

inline constexpr auto kReachable = reachability();

// The BGL model: a stateless tag, vertices = numbers, edges = table ids, the
// record as edge bundle (boost::get(&SubgroupRelation::index, g) is a property
// map). Parallel edges (distinct classes) and self-loops (isomorphic) occur.
struct TargetOf {
  [[nodiscard]] int operator()(int edge) const noexcept {
    return sub_of(relation(edge));
  }
};

struct SubgroupCsr {
  using vertex_descriptor = int;
  using edge_descriptor = int;
  using directed_category = boost::bidirectional_tag;
  using edge_parallel_category = boost::allow_parallel_edge_tag;
  struct traversal_category : boost::bidirectional_graph_tag,
                              boost::adjacency_graph_tag,
                              boost::vertex_list_graph_tag,
                              boost::edge_list_graph_tag {};
  using out_edge_iterator = int const *;
  using in_edge_iterator = int const *;
  using adjacency_iterator = boost::transform_iterator<TargetOf, int const *>;
  using vertex_iterator = boost::counting_iterator<int>;
  using edge_iterator = boost::counting_iterator<int>;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using degree_size_type = std::size_t;
  [[nodiscard]] static constexpr int null_vertex() noexcept { return -1; }

  using graph_bundled = boost::no_property;
  using vertex_bundled = boost::no_property;
  using edge_bundled = data::SubgroupRelation;
  [[nodiscard]] constexpr data::SubgroupRelation const &
  operator[](int edge) const noexcept {
    return relation(edge);
  }
};

using EdgeRange = std::pair<int const *, int const *>;

[[nodiscard]] constexpr EdgeRange as_range(std::span<int const> s) noexcept {
  return {s.data(), s.data() + s.size()};
}

[[nodiscard]] constexpr EdgeRange out_edges(int v,
                                            SubgroupCsr const &) noexcept {
  return as_range(kOutEdges[v]);
}
[[nodiscard]] constexpr std::size_t out_degree(int v,
                                               SubgroupCsr const &) noexcept {
  return kOutEdges[v].size();
}
[[nodiscard]] constexpr EdgeRange in_edges(int v,
                                           SubgroupCsr const &) noexcept {
  return as_range(kInEdges[v]);
}
[[nodiscard]] constexpr std::size_t in_degree(int v,
                                              SubgroupCsr const &) noexcept {
  return kInEdges[v].size();
}
[[nodiscard]] constexpr std::size_t degree(int v,
                                           SubgroupCsr const &g) noexcept {
  return out_degree(v, g) + in_degree(v, g);
}
[[nodiscard]] constexpr int source(int edge, SubgroupCsr const &) noexcept {
  return relation(edge).super;
}
[[nodiscard]] constexpr int target(int edge, SubgroupCsr const &) noexcept {
  return sub_of(relation(edge));
}
[[nodiscard]] inline std::pair<SubgroupCsr::adjacency_iterator,
                               SubgroupCsr::adjacency_iterator>
adjacent_vertices(int v, SubgroupCsr const &g) {
  auto const [first, last] = out_edges(v, g);
  return {SubgroupCsr::adjacency_iterator(first, TargetOf{}),
          SubgroupCsr::adjacency_iterator(last, TargetOf{})};
}
[[nodiscard]] inline std::pair<boost::counting_iterator<int>,
                               boost::counting_iterator<int>>
vertices(SubgroupCsr const &) noexcept {
  return {boost::counting_iterator<int>(0),
          boost::counting_iterator<int>(static_cast<int>(kNumVertices))};
}
[[nodiscard]] constexpr std::size_t num_vertices(SubgroupCsr const &) noexcept {
  return kNumVertices;
}
[[nodiscard]] inline std::pair<boost::counting_iterator<int>,
                               boost::counting_iterator<int>>
edges(SubgroupCsr const &) noexcept {
  return {boost::counting_iterator<int>(0),
          boost::counting_iterator<int>(static_cast<int>(kNumEdges))};
}
[[nodiscard]] constexpr std::size_t num_edges(SubgroupCsr const &) noexcept {
  return kNumEdges;
}

// Vertices are their own indices.
[[nodiscard]] inline boost::identity_property_map
get(boost::vertex_index_t, SubgroupCsr const &) noexcept {
  return {};
}

// The edge bundle's members as readable property maps over edge ids.
template <class T> struct Member {
  T data::SubgroupRelation::*member;
  [[nodiscard]] T const &operator()(int edge) const noexcept {
    return relation(edge).*member;
  }
};
template <class T>
using MemberMap = boost::function_property_map<Member<T>, int, T const &>;

template <class T>
[[nodiscard]] MemberMap<T> get(T data::SubgroupRelation::*member,
                               SubgroupCsr const &) noexcept {
  return MemberMap<T>(Member<T>{member});
}

// Edge predicate for the kind-restricted views: nullopt keeps every edge.
struct KindIs {
  std::optional<SubgroupKind> kind;
  [[nodiscard]] bool operator()(int edge) const noexcept {
    return !kind || relation(edge).kind == *kind;
  }
};

BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<SubgroupCsr>));
BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<SubgroupCsr>));
BOOST_CONCEPT_ASSERT((boost::EdgeListGraphConcept<SubgroupCsr>));
BOOST_CONCEPT_ASSERT((boost::AdjacencyGraphConcept<SubgroupCsr>));

} // namespace detail

// The one graph, and its kind-restricted view type. Both are ordinary BGL
// graphs: breadth_first_search, filtered_graph, reverse_graph and the visitors
// all apply to them unchanged.
inline constexpr detail::SubgroupCsr kSubgroupGraph{};
using SubgroupView = boost::filtered_graph<detail::SubgroupCsr, detail::KindIs>;

// The maximal-subgroup graph of the 230 space groups (t and k), constant-
// initialised from the table tools/generate_subgroup_relations.py derives.
class SubgroupGraph {
public:
  // Out-edge ids of `number`, of one kind or both; empty out of range.
  [[nodiscard]] static constexpr auto
  maximal_subgroups(int number,
                    std::optional<SubgroupKind> kind = std::nullopt) noexcept {
    return detail::kOutEdges[number] | std::views::filter(detail::KindIs{kind});
  }

  // In-edge ids of `number`.
  [[nodiscard]] static constexpr auto minimal_supergroups(
      int number, std::optional<SubgroupKind> kind = std::nullopt) noexcept {
    return detail::kInEdges[number] | std::views::filter(detail::KindIs{kind});
  }

  [[nodiscard]] static SubgroupEdge edge(int id) noexcept {
    return SubgroupEdge::of(id);
  }
  [[nodiscard]] static constexpr std::size_t num_edges() noexcept {
    return detail::kNumEdges;
  }

  // Reachability (reflexive), a bit test in the compile-time closure.
  [[nodiscard]] static constexpr bool is_subgroup(int sub, int super) noexcept {
    if (super < 1 || super > kNumSpaceGroups || sub < 1 ||
        sub > kNumSpaceGroups) {
      return false;
    }
    return sub == super ||
           detail::kReachable[static_cast<std::size_t>(super)].test(
               static_cast<std::size_t>(sub));
  }

  // A shortest edge chain super -> sub through steps of `kind` (both when
  // unset); empty for super == sub, nullopt when unreachable.
  [[nodiscard]] static std::optional<std::vector<int>>
  path(int super, int sub, std::optional<SubgroupKind> kind = std::nullopt);
};

} // namespace seitz::group

namespace boost {

template <>
struct property_map<seitz::group::detail::SubgroupCsr, vertex_index_t> {
  using type = identity_property_map;
  using const_type = identity_property_map;
};

template <class T>
struct property_map<seitz::group::detail::SubgroupCsr,
                    T seitz::data::SubgroupRelation::*> {
  using type = seitz::group::detail::MemberMap<T>;
  using const_type = type;
};

} // namespace boost

#pragma GCC visibility pop
