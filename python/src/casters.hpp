#pragma once

#include <cppcrystal/core/types.hpp>

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>

#include <cstddef>
#include <span>
#include <string>
#include <string_view>

// The conversions that are not pybind11's job out of the box: borrowing a span
// into an object's own storage, and the small NumPy shapes the library speaks.
namespace cppcrystal::python {

namespace py = pybind11;

// A span into `parent`'s own storage, as a Python list whose elements keep
// `parent` alive.
//
// Deliberately not a std::span type caster. A caster cannot see the parent, and
// copying the elements is the wrong answer for group::Wyckoff, whose ADDRESS is
// its identity -- generate::Placed stores a Wyckoff const *, so a copied
// Wyckoff would be an object whose address means nothing.
template <class T>
[[nodiscard]] py::list borrowed_list(py::handle parent,
                                     std::span<T const> items) {
  py::list out;
  for (T const &item : items) {
    out.append(
        py::cast(item, py::return_value_policy::reference_internal, parent));
  }
  return out;
}

// A span of values with no identity of their own, as a list of copies.
template <class T>
[[nodiscard]] py::list copied_list(std::span<T const> items) {
  py::list out;
  for (T const &item : items) {
    out.append(py::cast(item, py::return_value_policy::copy));
  }
  return out;
}

// std::string_view into a constexpr table, as a str. The copy is not an
// oversight to optimize away later: the tables outlive the interpreter, but a
// Python str has to own its bytes.
[[nodiscard]] inline py::str to_str(std::string_view text) {
  return py::str(text.data(), text.size());
}

// Cell::types() as an (N,) int32 array rather than a list[int], so it pairs
// with .positions for vectorized work instead of forcing a Python-level zip.
[[nodiscard]] inline py::array_t<int> types_array(Types const &types) {
  auto const count = static_cast<Index>(types.size());
  py::array_t<int> out(count);
  auto view = out.mutable_unchecked<1>();
  for (Index i = 0; i < count; ++i) {
    view(i) = types[static_cast<std::size_t>(i)];
  }
  return out;
}

} // namespace cppcrystal::python
