#pragma once

#include <boost/leaf.hpp>

#include <string>

namespace spglib {

namespace leaf = boost::leaf;
template <class T> using Result = boost::leaf::result<T>;

struct e_spacegroup_search_failed {};
struct e_cell_standardization_failed {};
struct e_symmetry_operation_search_failed {};
struct e_pointgroup_not_found {};
struct e_niggli_failed {};
struct e_delaunay_failed {};

struct e_atoms_too_close {
  double distance;
};

struct e_array_size_shortage {
  int provided;
  int needed;
};

// Free-form human-readable context, attached alongside a tag for diagnostics.
struct e_message {
  std::string text;
};

} // namespace spglib
