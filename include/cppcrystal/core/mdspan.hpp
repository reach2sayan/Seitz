#pragma once

// The one place the library says "mdspan". It binds to the vendored reference
// implementation (third_party/mdspan.hpp, Apache-2.0 WITH LLVM-exception,
// pinned in tools/mdspan_vendor_diff.py): the standard libraries this project
// builds against do not ship a complete <mdspan>, and a std::mdspan without
// std::submdspan leaves slicing with no spelling at all. Explicit `using`
// declarations, not a namespace alias: `namespace md = std` would make
// cppcrystal::md::vector compile.

// Before the include: otherwise the vendored header plants itself in ::std,
// which is UB and collides with a real <mdspan> in the same TU.
#ifndef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE cppcrystal_mdspan
#endif

// The vendored header is upstream-verbatim and not held to our warning set.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#include <cppcrystal/third_party/mdspan.hpp>
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include <cppcrystal/core/types.hpp>

#include <cstddef>

#pragma GCC visibility push(default)

namespace cppcrystal::md {

using MDSPAN_IMPL_STANDARD_NAMESPACE::dextents;
using MDSPAN_IMPL_STANDARD_NAMESPACE::dynamic_extent;
using MDSPAN_IMPL_STANDARD_NAMESPACE::extents;
using MDSPAN_IMPL_STANDARD_NAMESPACE::mdspan;

using MDSPAN_IMPL_STANDARD_NAMESPACE::layout_left;
using MDSPAN_IMPL_STANDARD_NAMESPACE::layout_right;
using MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride;

using MDSPAN_IMPL_STANDARD_NAMESPACE::full_extent;
using MDSPAN_IMPL_STANDARD_NAMESPACE::strided_slice;
using MDSPAN_IMPL_STANDARD_NAMESPACE::submdspan;

// A read-only view of a flat constexpr array as a fixed-shape table: the
// replacement for nested std::array<std::array<...>> lookup tables, usable in
// constant evaluation.
template <class T, std::size_t... Extents>
using table = mdspan<T const, extents<std::size_t, Extents...>>;

// A rows x cols view over a flat buffer (row-major), e.g. one permutation
// per symmetry operation.
template <class T> using matrix_view = mdspan<T, dextents<Index, 2>>;

} // namespace cppcrystal::md

#pragma GCC visibility pop
