#pragma once

#include <boost/leaf.hpp>

#include <concepts>
#include <string>
#include <type_traits>

#pragma GCC visibility push(default)

namespace cppcrystal {

namespace leaf = boost::leaf;
template <class T> using Result = boost::leaf::result<T>;

// A nullary callable that produces a Result. What the LEAF handling scopes
// take instead of a finished Result -- taking the call rather than its value
// is what lets a BOOST_LEAF_AUTO chain be written inline at the call site.
template <class F>
concept ResultProducer =
    std::invocable<F &> &&
    requires { typename std::invoke_result_t<F &>::value_type; };

struct e_spacegroup_search_failed {};
struct e_cell_standardization_failed {};
struct e_symmetry_operation_search_failed {};
struct e_magnetic_symmetry_search_failed {};
struct e_pointgroup_not_found {};
struct e_niggli_failed {};
struct e_delaunay_failed {};
struct e_empty_cell {};

// Degenerate input rejected at a public entry point: a (near-)singular lattice
// whose inverse would propagate NaN/inf through the whole pipeline. (A
// non-positive sampling mesh is unrepresentable rather than an error:
// kpoint::Mesh::of rejects it.)
struct e_invalid_lattice {
  double determinant;
};

struct e_atoms_too_close {
  double distance;
};

// Free-form human-readable context, attached alongside a tag for diagnostics.
struct e_message {
  std::string text;
};

} // namespace cppcrystal

#pragma GCC visibility pop
