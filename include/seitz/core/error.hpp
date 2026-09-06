#pragma once

#include <boost/leaf.hpp>

#include <concepts>
#include <cstdint>
#include <string>
#include <type_traits>

#pragma GCC visibility push(default)

namespace seitz {

namespace leaf = boost::leaf;
template <class T> using Result = boost::leaf::result<T>;

// A nullary callable producing a Result: what the LEAF handling scopes take
// instead of a finished Result, so a BOOST_LEAF_AUTO chain can be written
// inline at the call site.
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
// whose inverse would propagate NaN/inf. (A non-positive mesh is
// unrepresentable rather than an error -- kpoint::Mesh::of rejects it.)
struct e_invalid_lattice {
  double determinant;
};

struct e_atoms_too_close {
  double distance;
};

// A change of basis a Cell cannot take: a singular integer matrix, or one that
// mixes an aperiodic axis with a periodic one.
struct e_invalid_transformation {
  int determinant;
};

// A caller-supplied lattice whose metric the requested group's operations do
// not preserve, so no structure on it could carry that symmetry.
struct e_incompatible_lattice {};

// A CIF document the grammar could not consume to the end: the 1-based line
// and column the parse stopped at.
struct e_cif_syntax {
  std::int64_t line;
  std::int64_t column;
};

// A CIF tag the reader needs and the block does not carry.
struct e_cif_missing {
  std::string tag;
};

// A coordinate triplet ("x,-y,z+1/2") that is not three coordinates, or whose
// rotation part is not unimodular.
struct e_invalid_xyz {
  std::string text;
};

// A chemical symbol no tabulated element matches.
struct e_unknown_element {
  std::string symbol;
};

// A Hermann-Mauguin or Hall symbol no tabulated setting matches.
struct e_unknown_spacegroup_symbol {
  std::string symbol;
};

// Free-form human-readable context, attached alongside a tag for diagnostics.
struct e_message {
  std::string text;
};

} // namespace seitz

#pragma GCC visibility pop
