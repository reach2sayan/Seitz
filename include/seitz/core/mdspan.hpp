#pragma once

// The one place the library says "mdspan", bound to the vendored reference
// implementation (third_party/mdspan.hpp, Apache-2.0 WITH LLVM-exception,
// pinned in tools/mdspan_vendor_diff.py): no standard library built against
// here ships a complete <mdspan>, and mdspan without submdspan leaves slicing
// unspellable. Explicit `using`s, not a namespace alias -- `namespace md = std`
// would make seitz::md::vector compile.

// Before the include: otherwise the vendored header plants itself in ::std,
// which is UB and collides with a real <mdspan> in the same TU.
#ifndef MDSPAN_IMPL_STANDARD_NAMESPACE
#define MDSPAN_IMPL_STANDARD_NAMESPACE seitz_mdspan
#endif

// The vendored header is upstream-verbatim and not held to our warning set.
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#endif
#include <seitz/third_party/mdspan.hpp>
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#include <seitz/core/types.hpp>

#include <cstddef>

#pragma GCC visibility push(default)

namespace seitz::md {

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

// A flat constexpr array as a fixed-shape read-only table, replacing nested
// std::array lookup tables and usable in constant evaluation.
template <class T, std::size_t... Extents>
using table = mdspan<T const, extents<std::size_t, Extents...>>;

// A rows x cols view over a flat buffer (row-major), e.g. one permutation
// per symmetry operation.
template <class T> using matrix_view = mdspan<T, dextents<Index, 2>>;

} // namespace seitz::md

#pragma GCC visibility pop
